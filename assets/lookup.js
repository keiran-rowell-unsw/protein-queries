/* Preliminary browser-side UniProt identity lookup and immediate sequence summary.
 *
 * Retrieval: the browser (not the PUN) contacts rest.uniprot.org once per accession;
 * the JSON response already contains the canonical sequence, so no second fetch is
 * needed. Calculations: the PUN's /protein-prop-summary endpoint runs bounded,
 * offline, deterministic library calculations (no retrieval) on the supplied
 * sequence. Everything is rendered with text-safe DOM APIs; fetched JSON is never
 * inserted as HTML. Both displays are preliminary: the queued PBS job remains the
 * authoritative, verified report.
 */
(function () {
  "use strict";

  var UNIPROT = "https://rest.uniprot.org/uniprotkb/";
  var ACCESSION = /^([A-Za-z0-9]{6}|[A-Za-z0-9]{10})(-[0-9]{1,2})?$/;
  var input = document.getElementById("accession");
  var output = document.getElementById("lookup");
  var sequence = document.getElementById("sequence");
  var summarySection = document.getElementById("summary");
  var summaryButton = document.getElementById("summarize");
  var csrfField = document.getElementById("csrf");
  var prefix = document.body.getAttribute("data-prefix") || "";
  if (!input || !output || !sequence || !summarySection || !summaryButton || !csrfField) return;

  var timer = null;
  var request = 0;
  var canonical = ""; // canonical sequence from the once-off UniProt fetch

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
    paragraph(output, "Use “Immediate preliminary summary” for sequence properties, or “Queue full report” for the verified PBS report.");
  }

  function message(text) {
    clear(output);
    paragraph(output, text, "hint");
  }

  function showSummary(payload) {
    var props = payload.properties || {};
    clear(summarySection);
    summarySection.hidden = false;
    var heading = document.createElement("h2");
    heading.textContent = "Preliminary sequence summary";
    summarySection.appendChild(heading);
    paragraph(summarySection,
      "Computed now from the sequence alone; not retrieved, not verified against UniProt, "
      + "and not a substitute for the queued report.",
      "hint");
    if (identity.name) line(summarySection, "Protein", identity.name);
    if (identity.organism) line(summarySection, "Organism", identity.organism);
    if (props.length_residues) line(summarySection, "Length", props.length_residues + " residues");
    if (props.molecular_weight_da) line(summarySection, "Molecular weight", props.molecular_weight_da + " Da");
    if (props.isoelectric_point !== undefined) line(summarySection, "Isoelectric point", props.isoelectric_point);
    if (props.net_charge_pH_7 !== undefined) line(summarySection, "Net charge at pH 7", props.net_charge_pH_7);
    if (props.gravy_kyte_doolittle !== undefined) {
      line(summarySection, "Average hydropathy (GRAVY)", props.gravy_kyte_doolittle);
    }
    if (props.instability_index_guruprasad !== undefined) {
      line(summarySection, "Instability index", props.instability_index_guruprasad);
    }
    if (props["molar_extinction_coefficient_M-1cm-1"]) {
      var ext = props["molar_extinction_coefficient_M-1cm-1"];
      line(summarySection, "A280 (reduced / oxidised)", ext.reduced + " / " + ext.oxidised + " M⁻¹cm⁻¹");
    }
    if (props.inclusion_body_propensity) {
      var body = props.inclusion_body_propensity;
      line(summarySection, "Inclusion-body flag", String(body["class"] || "").replace("likely_", "likely "));
    }
    var warnings = props.warnings || [];
    for (var i = 0; i < warnings.length; i++) paragraph(summarySection, "Note: " + warnings[i], "hint");
    paragraph(summarySection,
      "The queued full report adds verified UniProt annotations, PDBe structures, "
      + "disorder prediction, and structural-neighbour search.");
  }

  var identity = {};

  function summarize() {
    var raw = sequence.value.trim();
    var payloadSequence = raw || canonical;
    if (!payloadSequence) {
      clear(summarySection);
      summarySection.hidden = false;
      paragraph(summarySection,
        "Enter a sequence, or look up an accession first, then run the immediate summary.", "hint");
      return;
    }
    var body = new URLSearchParams();
    body.set("csrf", csrfField.value);
    body.set("sequence", payloadSequence);
    summaryButton.disabled = true;
    fetch(prefix + "/protein-prop-summary", {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded" },
      body: body.toString(),
    })
      .then(function (response) {
        return response.json().then(function (data) {
          summaryButton.disabled = false;
          if (!response.ok) {
            clear(summarySection);
            summarySection.hidden = false;
            paragraph(summarySection, data.error || "The immediate summary failed.", "hint");
            return;
          }
          showSummary(data);
        });
      })
      .catch(function () {
        summaryButton.disabled = false;
        clear(summarySection);
        summarySection.hidden = false;
        paragraph(summarySection, "The immediate summary is unavailable; queue the full report.", "hint");
      });
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
          canonical = (entry.sequence && entry.sequence.value) || "";
          identity = { name: name, organism: organism };
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
  summaryButton.addEventListener("click", summarize);
})();
