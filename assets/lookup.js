/* Browser-side UniProt identity preview.
 *
 * Retrieval: the browser (not the PUN) contacts rest.uniprot.org once per typed
 * accession, so the submitter can confirm the right protein before queueing.
 * Everything is rendered with text-safe DOM APIs; fetched JSON is never inserted
 * as HTML. The preview is preliminary: the queued PBS job verifies identity (exact
 * sequence when one is pasted, accession precedence otherwise) and produces the
 * authoritative report and all of its figures on the compute node.
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

  function paragraph(node, text, className) {
    var note = document.createElement("p");
    if (className) note.className = className;
    note.textContent = text;
    node.appendChild(note);
    return note;
  }

  function showIdentity(preliminary) {
    clear(output);
    paragraph(output, "Preliminary browser lookup; the queued report is the verified result.", "hint");
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
    paragraph(output, "Queue the full report for verified annotations, structures, and figures.");
  }

  function message(text) {
    clear(output);
    paragraph(output, text, "hint");
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
          showIdentity({
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
