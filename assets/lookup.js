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
  var UNIPARC = "https://rest.uniprot.org/uniparc/search";
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
  var matchRequest = 0;
  var canonical = ""; // canonical sequence from the once-off UniProt fetch
  var matchSection = summarySection; // match result renders inside the summary card
  var matchStatus = document.createElement("div");

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

  function curve(container, title, points, xLabel, options) {
    var width = 640, height = 160, pad = 34;
    var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("viewBox", "0 0 " + width + " " + height);
    svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", title);
    svg.style.width = "100%";
    svg.style.height = "auto";
    var xs = points.map(function (p) { return p[0]; });
    var ys = points.map(function (p) { return p[1]; });
    var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
    var yMin = options.yMin !== undefined ? options.yMin : Math.min.apply(null, ys);
    var yMax = options.yMax !== undefined ? options.yMax : Math.max.apply(null, ys);
    if (yMax === yMin) yMax = yMin + 1;
    var px = function (x) { return pad + (x - xMin) / (xMax - xMin) * (width - pad - 10); };
    var py = function (y) { return height - pad - (y - yMin) / (yMax - yMin) * (height - pad - 14); };
    var polyline = document.createElementNS("http://www.w3.org/2000/svg", "polyline");
    polyline.setAttribute("fill", "none");
    polyline.setAttribute("stroke", options.color || "#0e7490");
    polyline.setAttribute("stroke-width", "2");
    polyline.setAttribute("points", points.map(function (p) {
      return px(p[0]).toFixed(1) + "," + py(p[1]).toFixed(1);
    }).join(" "));
    svg.appendChild(polyline);
    if (options.zero !== undefined) {
      var zero = document.createElementNS("http://www.w3.org/2000/svg", "line");
      zero.setAttribute("x1", pad); zero.setAttribute("x2", width - 10);
      zero.setAttribute("y1", py(options.zero)); zero.setAttribute("y2", py(options.zero));
      zero.setAttribute("stroke", "#94a3b8"); zero.setAttribute("stroke-dasharray", "4 4");
      svg.appendChild(zero);
      var zeroLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
      zeroLabel.setAttribute("x", 4); zeroLabel.setAttribute("y", py(options.zero) + 4);
      zeroLabel.setAttribute("font-size", "11"); zeroLabel.setAttribute("fill", "#64748b");
      zeroLabel.textContent = options.zeroLabel || "0";
      svg.appendChild(zeroLabel);
    }
    var caption = document.createElementNS("http://www.w3.org/2000/svg", "text");
    caption.setAttribute("x", pad); caption.setAttribute("y", height - 8);
    caption.setAttribute("font-size", "11"); caption.setAttribute("fill", "#64748b");
    caption.textContent = xLabel;
    svg.appendChild(caption);
    var heading = document.createElement("h3");
    heading.textContent = title;
    container.appendChild(heading);
    container.appendChild(svg);
    paragraph(container, options.note, "hint");
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
    if (payload.plot && payload.plot.points && payload.plot.points.length > 1) {
      curve(summarySection, "Hydropathy along the sequence", payload.plot.points,
        "Residue 1 – " + props.length_residues + " (" + payload.plot.scale + ", window "
        + payload.plot.window + ")",
        { zero: 0, zeroLabel: "0", color: "#0e7490",
          note: "Positive values are more hydrophobic; this is a plot of calculated "
            + "values, not a new measurement." });
    }
    var charge = props.charge_vs_ph && props.charge_vs_ph.ph_to_net_charge;
    if (charge) {
      var chargePoints = Object.keys(charge).map(function (pH) {
        return [parseFloat(pH), charge[pH]];
      }).sort(function (a, b) { return a[0] - b[0]; });
      curve(summarySection, "Net charge versus pH", chargePoints, "pH 3 – 11",
        { color: "#7c3aed", zero: 0, zeroLabel: "0",
          note: "Henderson–Hasselbalch model; charges near the pI are the least reliable." });
    }
    clear(matchStatus);
    matchSection.appendChild(matchStatus);
  }

  var identity = {};

  function showMatch(match) {
    clear(match.output);
    match.section.hidden = false;
    if (match.status === "verified") {
      var link = document.createElement("a");
      link.href = match.uniprot_url;
      link.rel = "noopener";
      link.textContent = match.accession;
      var heading = document.createElement("strong");
      heading.appendChild(document.createTextNode("Exact UniParc match: "));
      heading.appendChild(link);
      match.output.appendChild(heading);
      if (match.name) line(match.output, "Protein", match.name);
      if (match.organism) line(match.output, "Organism", match.organism);
      if (match.reviewed) line(match.output, "Reviewed", match.reviewed);
      paragraph(match.output,
        "The pasted sequence is identical to this UniProtKB entry (CRC64 "
        + match.checksum + ", verified residue by residue).", "hint");
    } else if (match.status === "ambiguous") {
      paragraph(match.output,
        "Several UniProtKB entries share this exact sequence. Enter the accession "
        + "you want to annotate. Candidates: " + match.candidates.join(", "), "hint");
    } else {
      paragraph(match.output, match.reason, "hint");
    }
  }

  function fetchEntry(accession, checksum, match) {
    return fetch(UNIPROT + encodeURIComponent(accession) + ".json")
      .then(function (response) {
        if (!response.ok) return null;
        return response.json();
      })
      .then(function (entry) {
        var description = (entry && entry.proteinDescription) || {};
        var recommended = description.recommendedName || {};
        match.accession = accession;
        match.checksum = checksum;
        match.uniprot_url = "https://www.uniprot.org/uniprotkb/" + accession + "/entry";
        match.name = (recommended.fullName || {}).value || "";
        match.organism = (((entry || {}).organism || {}).scientificName) || "";
        match.reviewed = entry && entry.entryType ? String(entry.entryType) : "";
        showMatch(match);
      });
  }

  function sequenceMatch(sequence, crc) {
    matchStatus.textContent = "Searching UniParc for an exact sequence match…";
    var current = ++matchRequest;
    fetch(
      UNIPARC + "?query=checksum%3A" + encodeURIComponent(crc)
      + "&format=json&size=5"
    )
      .then(function (response) {
        if (current !== matchRequest) return;
        if (!response.ok) {
          showMatch({ output: matchStatus, section: matchSection, status: "unavailable",
            reason: "The exact-sequence search is unavailable right now." });
          return;
        }
        var total = parseInt(response.headers.get("X-Total-Results") || "0", 10);
        return response.json().then(function (payload) {
          if (current !== matchRequest) return;
          if (total > 5) {
            showMatch({ output: matchStatus, section: matchSection, status: "ambiguous",
              candidates: ["too many to list"],
              reason: "Too many checksum hits to verify safely." });
            return;
          }
          var candidates = [];
          var results = (payload && payload.results) || [];
          for (var i = 0; i < results.length; i++) {
            var seq = results[i].sequence || {};
            if (seq.value !== sequence || seq.length !== sequence.length) continue;
            var accessions = results[i].uniProtKBAccessions || [];
            for (var j = 0; j < accessions.length; j++) {
              var acc = String(accessions[j]).split(".", 1)[0];
              if (/^[A-Z0-9]{6,10}$/.test(acc) && candidates.indexOf(acc) < 0) {
                candidates.push(acc);
              }
            }
          }
          if (!candidates.length) {
            showMatch({ output: matchStatus, section: matchSection, status: "not_found",
              reason: "No exact UniProt sequence match was found; identity is unresolved. "
                + "This search is exact, not a similarity search." });
            return;
          }
          if (candidates.length > 4) {
            showMatch({ output: matchStatus, section: matchSection, status: "ambiguous",
              candidates: candidates.slice(0, 8),
              reason: "Several UniProtKB accessions share this exact sequence." });
            return;
          }
          if (candidates.length === 1) {
            fetchEntry(candidates[0], crc,
              { output: matchStatus, section: matchSection, status: "verified" });
            return;
          }
          showMatch({ output: matchStatus, section: matchSection, status: "ambiguous",
            candidates: candidates, reason: "" });
        });
      })
      .catch(function () {
        if (current !== matchRequest) return;
        showMatch({ output: matchStatus, section: matchSection, status: "unavailable",
          reason: "The browser could not reach UniParc; identity match is unavailable." });
      });
  }

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
          sequenceMatch(data.sequence, data.crc64);
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
