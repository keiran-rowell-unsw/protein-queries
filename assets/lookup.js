/* Preliminary browser-side UniProt identity lookup.
 * The browser (not the PUN) contacts rest.uniprot.org, which allows cross-origin
 * GET requests. Everything is rendered with text-safe DOM APIs; fetched JSON is
 * never inserted as HTML. This display is preliminary: the queued PBS job is the
 * authoritative, verified report.
 */
(function () {
  "use strict";

  var UNIPROT = "https://rest.uniprot.org/uniprotkb/";
  var ACCESSION = /^([A-Za-z0-9]{6}|[A-Za-z0-9]{10})(-[0-9]{1,2})?$/;
  var input = document.getElementById("accession");
  var output = document.getElementById("lookup");
  if (!input || !output) return;

  var timer = null;
  var request = 0;

  function clear(node) {
    while (node.firstChild) node.removeChild(node.firstChild);
  }

  function line(node, label, value) {
    var row = document.createElement("div");
    var name = document.createElement("strong");
    name.textContent = label + ": ";
    row.appendChild(name);
    row.appendChild(document.createTextNode(value));
    node.appendChild(row);
  }

  function show(preliminary) {
    clear(output);
    var note = document.createElement("p");
    note.className = "hint";
    note.textContent = "Preliminary browser lookup; the queued report is the verified result.";
    output.appendChild(note);
    if (preliminary.uniprot_url) {
      var link = document.createElement("a");
      link.href = preliminary.uniprot_url;
      link.rel = "noopener";
      link.textContent = preliminary.uniprot_url;
      output.appendChild(link);
    }
    if (preliminary.name) line(output, "Protein", preliminary.name);
    if (preliminary.organism) line(output, "Organism", preliminary.organism);
    if (preliminary.length) line(output, "Length", preliminary.length + " residues");
    if (preliminary.reviewed) line(output, "Reviewed", preliminary.reviewed);
    var queueNote = document.createElement("p");
    queueNote.textContent = "Press “Queue report” to compute the full verified report on PBS.";
    output.appendChild(queueNote);
  }

  function message(text) {
    clear(output);
    var note = document.createElement("p");
    note.className = "hint";
    note.textContent = text;
    output.appendChild(note);
  }

  function lookup(accession) {
    if (!ACCESSION.test(accession)) {
      if (accession) message("Enter a 6- or 10-character UniProt accession.");
      else clear(output);
      return;
    }
    var current = ++request;
    fetch(UNIPROT + encodeURIComponent(accession) + ".json")
      .then(function (response) {
        if (current !== request) return;
        if (response.status === 404) {
          message("No UniProt entry found for " + accession + ".");
          return;
        }
        if (!response.ok) {
          message("UniProt returned an error (" + response.status + ").");
          return;
        }
        return response.json().then(function (entry) {
          if (current !== request) return;
          var primary = entry && entry.primaryAccession;
          if (!primary || primary !== accession.split("-")[0].toUpperCase()) {
            message("UniProt did not return the requested accession.");
            return;
          }
          var description = (entry.proteinDescription || {});
          var recommended = description.recommendedName || {};
          var name = (recommended.fullName || {}).value
            || (((description.submissionNames || [])[0] || {}).fullName || {}).value
            || "";
          var organism = (((entry.organism || {}).scientificName) || "");
          var reviewed = entry.entryType ? String(entry.entryType) : "";
          show({
            uniprot_url: "https://www.uniprot.org/uniprotkb/" + primary + "/entry",
            name: name,
            organism: organism,
            length: entry.sequence && entry.sequence.length,
            reviewed: reviewed,
          });
        });
      })
      .catch(function () {
        if (current !== request) return;
        message("The browser could not reach UniProt. You can still queue the report.");
      });
  }

  input.addEventListener("input", function () {
    var value = input.value.trim().toUpperCase();
    if (timer) window.clearTimeout(timer);
    timer = window.setTimeout(function () { lookup(value); }, 400);
  });
})();
