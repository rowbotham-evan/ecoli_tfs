# ecoli_tfs

A database of *Escherichia coli* transcription factors, published with GitHub Pages.

## Site

The site is plain static HTML/CSS/JS — no build step.

```
index.html         page layout
assets/style.css   styles (light + dark)
assets/app.js      loads data/dataset.csv into a searchable, sortable table
data/              put dataset.csv here
```

Preview locally:

```sh
python3 -m http.server 8000
# open http://localhost:8000
```

## Publishing

1. Push to `main` on GitHub.
2. In the repo: **Settings → Pages → Build and deployment → Deploy from a branch**, branch `main`, folder `/ (root)`.
3. The site will be live at `https://rowbotham-evan.github.io/ecoli_tfs/`.
