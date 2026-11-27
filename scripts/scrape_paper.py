#!/usr/bin/env python3
"""Scrape one open-access paper and catalogue it in data/ingestion/papers.json.

    python3 scripts/scrape_paper.py DOI --method PPTP
    python3 scripts/scrape_paper.py --batch papers.json     # a list of {"doi", "method"}

The paper is looked up by DOI in Europe PMC, and its full text is fetched as JATS XML
from the Europe PMC REST API (open-access articles only). From the XML the script takes
the bibliographic details, the section outline, the Methods section split into its
subsections, the figure, reference and supplementary-file counts, the data and code
accessions, the counts the paper states about its experiment (TFs, promoters, media), the
strains and growth conditions named in the Methods, and which of the database's TFs and
catalogued methods the text mentions. Each paper's entry also summarizes the records
data/dataset.csv already holds from it (TFs, target genes, media, response direction).
It also records how the scrape itself went (bytes, seconds, checksum), which is what
the Ingestion page shows.

The Methods text is stored only when the article's licence is CC BY or CC0, which allow
redistribution with attribution; otherwise only headings and word counts are kept.

--method is the dataset's method code the paper is the source for (SELEX, DAP, CHIP,
PPTP), or any other label for a method that has no records yet.

Scraping only writes the catalogue; it never adds records to data/dataset.csv. A paper
with no records there is marked "awaiting approval" until its data is imported.

Repeated runs replace the paper's entry (matched by DOI) and leave the others alone.
Uses only the standard library. The XML is read with html.parser, a pure-Python parser,
so the script does not depend on the interpreter's expat build.
"""

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATASET = ROOT / "data" / "dataset.csv"
CATALOGUE = ROOT / "data" / "ingestion" / "papers.json"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest"
SCRAPER_VERSION = 1
USER_AGENT = "paper-scraper/1 (python-urllib)"

# The methods catalogued in notes.md, and how papers write them.
METHOD_PATTERNS = {
    "gSELEX": r"\b(?:g|genomic )?SELEX\b",
    "DAP-seq": r"\b(?:multi)?DAP-seq\b",
    "ChIP-seq": r"\bChIP-(?:seq|exo|chip)\b",
    "PPTP-seq": r"\bPPTP-seq\b",
    "IPOD-HR": r"\bIPOD(?:-HR)?\b",
}
ACCESSION_PATTERNS = [
    ("GEO", r"\bGSE\d{3,}\b", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}"),
    ("BioProject", r"\bPRJ[NED][A-Z]\d+\b", "https://www.ncbi.nlm.nih.gov/bioproject/{}"),
    ("SRA", r"\b[SED]R[PRX]\d{5,}\b", "https://www.ncbi.nlm.nih.gov/sra/{}"),
    ("GenBank", r"\b(?:U00096|NC_000913)\.\d\b", "https://www.ncbi.nlm.nih.gov/nuccore/{}"),
    ("Zenodo", r"\b10\.5281/zenodo\.\d+\b", "https://doi.org/{}"),
    ("GitHub", r"\bgithub\.com/[\w.-]+/[\w.-]+", "https://{}"),
]
# Top-level sections that are not part of the article's main text.
BACK_TYPES = {"ack", "fn-group", "ref-list", "associated-data", "supplementary-material"}
BACK_TITLES = re.compile(
    r"^(acknowledge?ments?|author contributions?|peer review|competing interests?|footnotes|references|"
    r"supplementary (information|materials?|data)|associated data|funding( statement| information)?|abbreviations|"
    r"(data|code)( and code)? availability( statement)?|conflicts? of interests?|institutional review board statement|"
    r"informed consent statement|author summary|supporting information|declarations?|ethics statement)$", re.I)
METHODS_TITLE = re.compile(r"^(online |star[ \u2605]*)?(materials? (and|&) )?methods( and materials)?( details?)?$|"
                           r"^experimental( procedures?| section)?$|^methodology$", re.I)
# Counts a paper reports about its experiment: a number, up to three words, then one of these nouns.
NUMBER_WORDS = "one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve"
COUNTED = (r"TF genes|TFs|transcription factors|regulators|promoters|operons|sgRNAs?(?: sequences)?|media|strains|"
           r"binding sites|binding regions|peaks|uncharacterized proteins|DNA-binding proteins|candidate proteins|genomes|bins|biological replicates|replicates|conditions|"
           r"TF-gene responses|responses|experiments|motifs")
REPORTED = re.compile(rf"(?<![\w.])((?:more than|over|about|at least|nearly|up to)\s+)?((?:\d[\d,]*|{NUMBER_WORDS}))\s+((?:(?!and\b|or\b|of\b|in\b)[\w.-]+\s+){{0,3}}?(?:{COUNTED}))\b", re.I)
# Experimental conditions looked for in the Methods text.
CONDITION_PATTERNS = {
    "Strains": r"\b(?:MG1655|W3110|BW25113|DY330|FR-E01|DH5α|NEB 10-β|BL21(?:\(DE3\))?)(?![\w-])",
    "Growth media": r"\bM9 minimal (?:glucose|glycerol)(?: \([\d.]+%\))?|\bM9 (?:glucose|glycerol)\b|\bLB\b|\bMOPS\b|\bSOC\b",
    "Temperatures": r"\b\d{2}\s?°C",
    "Inducers and antibiotics": r"\b(?:aTc|anhydrotetracycline|IPTG|arabinose|kanamycin|ampicillin|chloramphenicol|spectinomycin|formaldehyde)\b",
    "Software": r"\b(?:Bowtie2?|BEDTools|SAMtools|MACS[23]?|MEME|Tomtom|DESeq2|Cutadapt|Trimmomatic|BWA|Scipy|FlowJo)\b(?: v?\d[\d.]*\d)?",
}
BLOCKS = {"p", "title", "sec", "label", "caption", "disp-formula", "list-item", "td", "th"}


class Node:
    def __init__(self, tag, attrs, parent=None):
        self.tag, self.attrs, self.parent, self.kids = tag, dict(attrs), parent, []

    def iter(self, tag=None):
        if tag is None or self.tag == tag:
            yield self
        for kid in self.kids:
            if isinstance(kid, Node):
                yield from kid.iter(tag)

    def children(self, tag):
        return [k for k in self.kids if isinstance(k, Node) and k.tag == tag]

    def child(self, tag):
        return next(iter(self.children(tag)), None)

    def text(self, skip=("fig", "table-wrap", "supplementary-material")):
        """The element's text. Citation markers are dropped, an equation becomes a
        placeholder, and block elements are kept apart by a space."""
        out = []

        def walk(node):
            for kid in node.kids:
                if isinstance(kid, str):
                    out.append(kid)
                elif kid.tag in skip:
                    continue
                elif kid.tag == "sup" and any(x.attrs.get("ref-type") == "bibr" for x in kid.iter("xref")):
                    continue
                elif kid.tag == "disp-formula":
                    out.append(" [equation] ")
                else:
                    walk(kid)
                    if kid.tag in BLOCKS:
                        out.append(" ")

        walk(self)
        return " ".join("".join(out).split())


class TreeBuilder(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root = self.cur = Node("root", {})

    def handle_starttag(self, tag, attrs):
        node = Node(tag, attrs, self.cur)
        self.cur.kids.append(node)
        self.cur = node

    def handle_startendtag(self, tag, attrs):
        self.cur.kids.append(Node(tag, attrs, self.cur))

    def handle_endtag(self, tag):
        node = self.cur
        while node is not self.root and node.tag != tag:
            node = node.parent
        if node is not self.root:
            self.cur = node.parent

    def handle_data(self, data):
        self.cur.kids.append(data)


def parse(xml):
    builder = TreeBuilder()
    builder.feed(xml)
    builder.close()
    return builder.root


def fetch(url, timeout=240, tries=4):
    """Returns (status, body bytes, seconds taken). A server error (Europe PMC answers 503
    when busy) or a timeout is retried after a pause."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(1, tries + 1):
        start = time.monotonic()
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                body = response.read()
                return response.status, body, round(time.monotonic() - start, 1)
        except urllib.error.HTTPError as error:
            if error.code < 500 or attempt == tries:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries:
                raise
        time.sleep(10 * attempt)


def words(text):
    return len(text.split())


def title_of(sec):
    title = sec.child("title")
    # Some journals number their headings ("4. Materials and Methods", "4.1. Strains").
    return re.sub(r"^\d+(\.\d+)*\.?\s+", "", title.text()) if title else ""


def section_text(sec):
    """A section's own paragraphs (not its title or subsections), one per line."""
    return "\n\n".join(p.text() for p in sec.children("p"))


def tf_names():
    """Distinct TF names in the database, with compound names such as YciT(DeoT) split."""
    with DATASET.open(newline="") as f:
        names = {r["tf_name"] for r in csv.DictReader(f)}
    return {name: [part for part in re.split(r"[()/\s]+", name) if len(part) >= 3] for name in names}


def count_tfs(text):
    """How often each database TF is named, as the protein (ArcA) or the gene (arcA)."""
    found = {}
    for name, parts in tf_names().items():
        n = sum(len(re.findall(rf"(?<![\w-])[{p[0].upper()}{p[0].lower()}]{re.escape(p[1:])}(?![\w])", text)) for p in parts)
        if n:
            found[name] = n
    # The methods the database already has records for, per TF (the base table's blank method aside).
    have = {}
    with DATASET.open(newline="") as f:
        for r in csv.DictReader(f):
            if r["method"]:
                have.setdefault(r["tf_name"], set()).add(r["method"])
    return [{"tf": tf, "mentions": n, "methods": sorted(have.get(tf, ()))}
            for tf, n in sorted(found.items(), key=lambda kv: (-kv[1], kv[0]))]


def read_tables(root, keep_cells):
    """The article's tables: label, caption, column names and size, and the cells when the
    licence allows storing them."""
    tables = []
    for wrap in root.iter("table-wrap"):
        rows = [[cell.text() for cell in tr.kids if isinstance(cell, Node) and cell.tag in ("td", "th")] for tr in wrap.iter("tr")]
        rows = [r for r in rows if any(r)]
        if not rows:
            continue
        label, caption = wrap.child("label"), wrap.child("caption")
        table = {"label": label.text().rstrip(".") if label else "", "caption": caption.text() if caption else "",
                 "columns": rows[0], "n_rows": len(rows) - 1}
        if keep_cells:
            table["rows"] = rows[1:]
        tables.append((table, rows[1:]))
    return tables


def studied_tfs(tables):
    """The TFs a paper studied, read from the table that lists the most gene names in its
    first column, each matched to the database: its TF name there and the methods that
    already have records for it."""
    gene = re.compile(r"^([a-z]{3}[A-Z]?)\b")
    best = max(tables, key=lambda t: sum(1 for r in t[1] if r and gene.match(r[0])), default=None)
    genes = [gene.match(r[0]).group(1) for r in best[1] if r and gene.match(r[0])] if best else []
    if len(genes) < 5:
        return None
    by_part = {part.lower(): name for name, parts in tf_names().items() for part in parts}
    have = {}
    with DATASET.open(newline="") as f:
        for r in csv.DictReader(f):
            if r["method"]:
                have.setdefault(r["tf_name"], set()).add(r["method"])
    tfs = []
    for g in dict.fromkeys(genes):
        name = by_part.get(g.lower(), "")
        tfs.append({"gene": g, "tf": name, "methods": sorted(have.get(name, ()))})
    return {"table": best[0]["label"], "caption": best[0]["caption"], "tfs": tfs}


def reported_numbers(parts):
    """Counts the paper states about its experiment. parts is [(where, text)]; each hit keeps
    its sentence so the number can be read in context."""
    found, seen = [], set()
    for where, text in parts:
        for sentence in re.split(r"(?<=[.!?])\s+(?=[A-Z])", text):
            for qualifier, value, what in REPORTED.findall(sentence):
                key = (value.lower(), what.lower())
                if key not in seen:
                    seen.add(key)
                    found.append({"qualifier": qualifier.strip().lower(), "value": value, "what": " ".join(what.split()), "where": where, "sentence": sentence.strip()})
    return found


def conditions(text):
    """Strains, media and so on named in the Methods, most mentioned first."""
    out = {}
    for group, pattern in CONDITION_PATTERNS.items():
        counts = {}
        for hit in re.findall(pattern, text):
            hit = hit.replace("\u2009", " ").replace(" °C", "°C")
            counts[hit] = counts.get(hit, 0) + 1
        if counts:
            out[group] = [{"term": t, "mentions": n} for t, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    return out


def database_summary(paper):
    """What data/dataset.csv holds from this paper: its records are those whose source
    names the paper's DOI or PMCID."""
    keys = [k.lower() for k in (paper.get("doi"), paper.get("pmcid")) if k]
    unsourced = paper.get("records_method")   # e.g. SELEX: the TEC records name no source
    with DATASET.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if any(k in r["source"].lower() for k in keys)
                or (unsourced and r["method"] == unsourced and not r["source"])]
    if not rows:
        return {"records": 0}
    per_tf, media, strains = {}, {}, {}
    for r in rows:
        per_tf[r["tf_name"]] = per_tf.get(r["tf_name"], 0) + 1
        direction = "up" if r["response"].startswith("up") else "down" if r["response"].startswith("down") else "other"
        m = media.setdefault(r["growth_medium"] or "Not recorded", {"records": 0, "up": 0, "down": 0, "other": 0})
        m["records"] += 1
        m[direction] += 1
        strains[r["strain"] or "Not recorded"] = strains.get(r["strain"] or "Not recorded", 0) + 1
    ranked = sorted(per_tf.items(), key=lambda kv: (-kv[1], kv[0]))
    return {
        "records": len(rows),
        "tfs": len(per_tf),
        "target_genes": len({r["associated_genes"] for r in rows if r["associated_genes"]}),
        "up": sum(m["up"] for m in media.values()),
        "down": sum(m["down"] for m in media.values()),
        "media": [{"medium": k, **v} for k, v in sorted(media.items(), key=lambda kv: -kv[1]["records"])],
        "strains": [{"strain": k, "records": n} for k, n in sorted(strains.items(), key=lambda kv: -kv[1])],
        "records_per_tf_median": sorted(per_tf.values())[len(per_tf) // 2],
        "top_tfs": [{"tf": tf, "records": n} for tf, n in ranked[:15]],
    }


def scrape(doi, method):
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi.strip())
    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # 1. Look the DOI up.
    query = urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json", "resultType": "core"})
    _, body, lookup_seconds = fetch(f"{API}/search?{query}")
    hits = json.loads(body)["resultList"]["result"]
    if not hits:
        raise SystemExit(f"Europe PMC has no record for DOI {doi}.")
    hit = hits[0]
    pmcid = hit.get("pmcid")
    if not pmcid or hit.get("isOpenAccess") != "Y":
        raise SystemExit(f"{doi} has no open-access full text in Europe PMC (PMCID: {pmcid or 'none'}).")
    licence = (hit.get("license") or "").lower()
    may_store_text = licence in ("cc by", "cc0")

    # 2. Fetch the full text.
    url = f"{API}/{pmcid}/fullTextXML"
    status, raw, seconds = fetch(url)
    xml = raw.decode("utf-8")
    root = parse(xml)
    article, body_el = next(root.iter("article")), next(root.iter("body"))
    meta = next(root.iter("article-meta"))

    # 3. Bibliographic details.
    authors = [f"{c.child('name').child('given-names').text()} {c.child('name').child('surname').text()}"
               for c in meta.iter("contrib") if c.child("name")]
    surnames = [a.split()[-1] for a in authors]
    pub = meta.child("pub-date")
    year = pub.child("year").text() if pub and pub.child("year") else hit.get("pubYear", "")
    lead = surnames[0] if len(surnames) == 1 else f"{surnames[0]} & {surnames[1]}" if len(surnames) == 2 else f"{surnames[0]} et al."
    abstract = next((a for a in meta.children("abstract") if "abstract-type" not in a.attrs), None)
    abstract_text = " ".join(p.text() for p in abstract.iter("p")) if abstract else ""

    # 4. Section outline. Europe PMC puts the back matter in the body too, so it is told apart here.
    outline, main_secs, methods_sec = [], [], None
    for sec in body_el.children("sec"):
        title = title_of(sec)
        if not title:
            continue
        is_back = sec.attrs.get("sec-type") in BACK_TYPES or bool(BACK_TITLES.match(title))
        is_methods = not is_back and bool(METHODS_TITLE.match(title))
        if is_methods:
            methods_sec = sec
        if not is_back:
            main_secs.append(sec)
            outline.append({"title": title, "words": words(sec.text()) - words(title),
                            "subsections": len(sec.children("sec")), "kind": "methods" if is_methods else "text"})
    main_text = " ".join(sec.text() for sec in main_secs)
    captions = " ".join(fig.text(skip=()) for fig in root.iter("fig"))

    # 5. The Methods section, one entry per subsection.
    methods, methods_text = [], []
    if methods_sec:
        subsections = methods_sec.children("sec") or [methods_sec]
        for sub in subsections:
            text = "\n\n".join(filter(None, [section_text(sub)] + [f"{title_of(s)}. {section_text(s)}" for s in sub.children("sec")]))
            methods_text.append(text)
            entry = {"heading": title_of(sub), "words": words(text), "paragraphs": text.count("\n\n") + 1 if text else 0,
                     "equations": len(list(sub.iter("disp-formula")))}
            if may_store_text:
                entry["text"] = text
            methods.append(entry)

    for entry in outline:   # the same count as the subsections' total, which leaves out their headings
        if entry["kind"] == "methods":
            entry["words"] = sum(m["words"] for m in methods)

    # 6. Supplementary files (Europe PMC lists them twice: where cited and under Associated Data).
    supplements, seen = [], set()
    for supp in root.iter("supplementary-material"):
        media = next(supp.iter("media"), None)
        name = media.attrs.get("xlink:href", "") if media else ""
        if name and name not in seen:
            seen.add(name)
            supplements.append({"label": supp.text(skip=()) or name, "file": name, "type": name.rsplit(".", 1)[-1].lower()})

    # 7. Data and code accessions, from the whole article outside the reference list.
    availability = " ".join(sec.text() for sec in body_el.children("sec") if re.search(r"availability", title_of(sec), re.I))
    searchable = " ".join([main_text, availability, captions])
    accessions = []
    for kind, pattern, link in ACCESSION_PATTERNS:
        for acc in sorted(set(re.findall(pattern, searchable))):
            acc = acc.rstrip(".")
            accessions.append({"kind": kind, "id": acc, "url": link.format(acc), "mentions": searchable.count(acc)})

    tables = read_tables(root, may_store_text)
    refs = list(root.iter("ref"))
    tfs = count_tfs(" ".join([main_text, captions]))
    method_mentions = {name: len(re.findall(pattern, searchable)) for name, pattern in METHOD_PATTERNS.items()}
    methods_words = sum(m["words"] for m in methods)

    return {
        "id": pmcid,
        "status": "scraped",
        "method": method,
        "citation": f"{lead} {year}",
        "title": meta.child("title-group").child("article-title").text(),
        "authors": authors,
        "journal": next(root.iter("journal-title")).text(),
        "year": int(year) if str(year).isdigit() else year,
        "volume": meta.child("volume").text() if meta.child("volume") else "",
        "article_type": article.attrs.get("article-type", ""),
        "doi": doi,
        "pmid": hit.get("pmid", ""),
        "pmcid": pmcid,
        "licence": licence.upper() or "unknown",
        "open_access": True,
        "scrape": {
            "source": "Europe PMC full-text XML (JATS)",
            "url": url,
            "retrieved": retrieved,
            "http_status": status,
            "bytes": len(raw),
            "seconds": seconds,
            "lookup_seconds": lookup_seconds,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "scraper": f"scripts/scrape_paper.py v{SCRAPER_VERSION}",
            "methods_text_stored": may_store_text,
        },
        "stats": {
            "xml_elements": sum(1 for _ in root.iter()) - 1,
            "main_sections": len(main_secs),
            "main_words": words(main_text),
            "abstract_words": words(abstract_text),
            "methods_subsections": len(methods),
            "methods_words": methods_words,
            "methods_share": round(methods_words / words(main_text), 3) if main_text else 0,
            "figures": len(list(root.iter("fig"))),
            "tables": len(list(root.iter("table-wrap"))),
            "equations": len(list(root.iter("disp-formula"))),
            "references": len(refs),
            "references_with_doi": sum(1 for r in refs if any(e.attrs.get("ext-link-type") == "doi" for e in r.iter("ext-link"))),
            "supplementary_files": len(supplements),
            "accessions": len(accessions),
            "tfs_mentioned": len(tfs),
            "tfs_in_database": len(tf_names()),
            "tf_mentions": sum(t["mentions"] for t in tfs),
        },
        # What the scraper looks for in every paper, and whether this one had it.
        "fields": {
            "Title and authors": bool(authors),
            "Abstract": bool(abstract_text),
            "Methods section": bool(methods),
            "Methods subsections": len(methods) > 1,
            "Data availability": bool(availability),
            "Data accessions": bool(accessions),
            "Supplementary files": bool(supplements),
            "Reference list": bool(refs),
            "Tables in the text": bool(list(root.iter("table-wrap"))),
        },
        "reported": reported_numbers([("Abstract", abstract_text)] + [(m["heading"], t) for m, t in zip(methods, methods_text)]),
        "conditions": conditions("\n".join(methods_text)),
        "studied": studied_tfs(tables),
        "tables": [t for t, _ in tables],
        "outline": outline,
        "methods": methods,
        "supplements": supplements,
        "accessions": accessions,
        "method_mentions": method_mentions,
        "tfs": tfs,
    }


def save(record, records_method=None):
    """Adds or replaces the paper's entry (matched by DOI), refreshes every entry's database
    figures, and writes the catalogue. Returns whether the paper was new."""
    CATALOGUE.parent.mkdir(parents=True, exist_ok=True)
    catalogue = json.loads(CATALOGUE.read_text()) if CATALOGUE.exists() else {"papers": []}
    papers = catalogue["papers"]
    at = next((i for i, p in enumerate(papers) if p.get("doi", "").lower() == record["doi"].lower()), None)
    records_method = records_method or (papers[at].get("records_method") if at is not None else None)
    if records_method:
        record["records_method"] = records_method
    if at is not None:   # the lab and the TFs characterized are curated by hand, so a rescrape keeps it
        for key in ("lab", "lab_group", "institution", "characterized"):
            if key in papers[at]:
                record.setdefault(key, papers[at][key])
    if at is None:
        papers.append(record)
    else:
        papers[at] = record
    for paper in papers:   # every paper's database figures, scraped or not, follow the current dataset
        paper["database"] = database_summary(paper)
        # Scraping never touches data/dataset.csv. A paper with no records there is a
        # candidate, shown on the Ingestion page until its data is approved and imported.
        paper["integration"] = "in database" if paper["database"]["records"] else "awaiting approval"
    CATALOGUE.write_text(json.dumps(catalogue, indent=1, ensure_ascii=False) + "\n")
    return at is None


def report(record, new):
    s, f, d = record["stats"], record["scrape"], record["database"]
    print(f"{record['citation']} ({record['pmcid']}): {f['bytes']:,} bytes in {f['seconds']} s")
    print(f"  {s['main_words']:,} words in {s['main_sections']} sections; Methods: {s['methods_words']:,} words in {s['methods_subsections']} subsections")
    print(f"  {s['figures']} figures, {s['references']} references, {s['supplementary_files']} supplementary files, {s['accessions']} accessions")
    print(f"  {s['tfs_mentioned']} of {s['tfs_in_database']} database TFs mentioned")
    print(f"  in the database: {d['records']:,} records" + (f", {d['tfs']} TFs, {d['target_genes']:,} target genes" if d["records"] else ""))
    print(f"  {'added to' if new else 'updated in'} {CATALOGUE.relative_to(ROOT)}")


def run_batch(path, workers=4):
    """Scrapes every paper in a JSON list of {"doi", "method"} objects, a few at a time.
    Papers already scraped are skipped; one that fails is reported and the rest carry on."""
    wanted = json.loads(Path(path).read_text())
    catalogue = json.loads(CATALOGUE.read_text()) if CATALOGUE.exists() else {"papers": []}
    done = {p["doi"].lower() for p in catalogue["papers"] if p.get("status") == "scraped"}
    todo, seen = [], set()
    for item in wanted:
        doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", item["doi"].strip()).lower()
        if doi not in done and doi not in seen:
            seen.add(doi)
            todo.append((doi, item["method"]))
    print(f"{len(wanted)} listed, {len(wanted) - len(todo)} already scraped or repeated, {len(todo)} to scrape", flush=True)
    failed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(scrape, doi, method): doi for doi, method in todo}
        for n, job in enumerate(as_completed(jobs), 1):
            doi = jobs[job]
            try:
                record = job.result()
            except BaseException as error:   # SystemExit included: no open full text
                failed.append((doi, str(error) or type(error).__name__))
                print(f"[{n}/{len(todo)}] FAILED {doi}: {failed[-1][1]}", flush=True)
                continue
            save(record)   # written from this thread only, one paper at a time
            st = record["stats"]
            print(f"[{n}/{len(todo)}] {record['citation']} ({record['method']}): Methods {st['methods_words']:,} words, {st['tfs_mentioned']} TFs named", flush=True)
    print(f"scraped {len(todo) - len(failed)}, failed {len(failed)}")
    for doi, why in failed:
        print(f"  {doi}: {why}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("doi", nargs="?")
    parser.add_argument("--method", help="method code the paper is the source for, e.g. PPTP")
    parser.add_argument("--records-method", help="also count the dataset's records of this method that name no source "
                        "(SELEX: the TEC database records) as this paper's")
    parser.add_argument("--batch", help='scrape every paper in this JSON file, a list of {"doi", "method"} objects')
    args = parser.parse_args()
    if args.batch:
        return run_batch(args.batch)
    if not args.doi or not args.method:
        parser.error("give a DOI and --method, or --batch FILE")
    record = scrape(args.doi, args.method)
    new = save(record, args.records_method)
    report(record, new)


if __name__ == "__main__":
    main()
