// Loads data/dataset.csv and renders a searchable, sortable table.
const DATA_URL = "data/dataset.csv";
const MAX_ROWS = 500;
const METHODS = ["SELEX", "TTPT", "DAP", "CHIP"];
// Dropdown edits live in this browser only; "Download CSV" exports them.
const EDITS_KEY = "method-edits-v2";

function parseCSV(text) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === ",") { row.push(field); field = ""; }
    else if (c === "\n" || c === "\r") {
      if (c === "\r" && text[i + 1] === "\n") i++;
      row.push(field); rows.push(row); row = []; field = "";
    } else field += c;
  }
  if (field || row.length) { row.push(field); rows.push(row); }
  return rows.filter(r => r.some(v => v !== ""));
}

const table = document.getElementById("table");
const search = document.getElementById("search");
const count = document.getElementById("count");
const empty = document.getElementById("empty");

let header = [], rows = [], sortCol = -1, sortDir = 1, methodCol = -1;

function loadEdits() {
  try { return JSON.parse(localStorage.getItem(EDITS_KEY)) || {}; } catch { return {}; }
}

function saveEdit(id, value) {
  const edits = loadEdits();
  edits[id] = value;
  try { localStorage.setItem(EDITS_KEY, JSON.stringify(edits)); } catch {}
}

function cell(r, v, i) {
  if (i !== methodCol) {
    return `<td title="${escape(v)}">${escape(v.length > 60 ? v.slice(0, 60) + "…" : v)}</td>`;
  }
  const options = ["", ...METHODS, ...(v && !METHODS.includes(v) ? [v] : [])];
  return `<td><select data-id="${escape(r[0])}">` + options.map(o =>
    `<option${o === v ? " selected" : ""}>${escape(o)}</option>`).join("") + "</select></td>";
}

function toCSV(data) {
  return data.map(r => r.map(v => /[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v).join(",")).join("\n") + "\n";
}

function render() {
  const q = search.value.trim().toLowerCase();
  let view = q ? rows.filter(r => r.some(v => v.toLowerCase().includes(q))) : rows.slice();
  if (sortCol >= 0) {
    view.sort((a, b) => sortDir * a[sortCol].localeCompare(b[sortCol], undefined, { numeric: true }));
  }
  table.tHead.innerHTML = "<tr>" + header.map((h, i) =>
    `<th data-i="${i}">${h}${i === sortCol ? (sortDir > 0 ? " ▲" : " ▼") : ""}</th>`).join("") + "</tr>";
  const body = view.slice(0, MAX_ROWS).map(r =>
    "<tr>" + r.map((v, i) => cell(r, v, i)).join("") + "</tr>");
  table.tBodies[0].innerHTML = body.join("");
  count.textContent = `${view.length} of ${rows.length} entries` + (view.length > MAX_ROWS ? ` (showing first ${MAX_ROWS})` : "");
}

function escape(s) {
  return s.replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

table.tHead.addEventListener("click", e => {
  const th = e.target.closest("th");
  if (!th) return;
  const i = Number(th.dataset.i);
  sortDir = i === sortCol ? -sortDir : 1;
  sortCol = i;
  render();
});
search.addEventListener("input", render);

table.tBodies[0].addEventListener("change", e => {
  const id = e.target.dataset.id;
  const row = rows.find(r => r[0] === id);
  if (!row) return;
  row[methodCol] = e.target.value;
  saveEdit(id, e.target.value);
});

document.getElementById("download").addEventListener("click", () => {
  const blob = new Blob([toCSV([header, ...rows])], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "dataset.csv";
  a.click();
  URL.revokeObjectURL(a.href);
});

fetch(DATA_URL)
  .then(r => { if (!r.ok) throw new Error(r.status); return r.text(); })
  .then(text => {
    [header = [], ...rows] = parseCSV(text);
    if (!header.length) throw new Error("empty");
    methodCol = header.indexOf("method");
    const edits = loadEdits();
    if (methodCol >= 0) rows.forEach(r => { if (r[0] in edits) r[methodCol] = edits[r[0]]; });
    render();
    empty.hidden = rows.length > 0;
  })
  .catch(() => {
    empty.hidden = false;
    table.hidden = true;
    search.disabled = true;
  });
