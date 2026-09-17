## report.R -------------------------------------------------------------------
## A dependency-free HTML index of everything a run produced. Enough to send a
## collaborator a single link instead of a folder of PNGs.
## -----------------------------------------------------------------------------

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' Build results/<contrast>/report.html linking every figure and table.
build_report <- function(cfg, cohort_summary = NULL, de_table = NULL, dmp_table = NULL) {
  figs <- sort(list.files(cfg$figures_dir, pattern = "\\.(png|jpg|svg)$", full.names = FALSE))
  tabs <- sort(list.files(cfg$tables_dir, pattern = "\\.tsv$", full.names = FALSE))

  fmt_table <- function(df, max_rows = 25) {
    if (is.null(df) || nrow(df) == 0) return("<p><em>no rows</em></p>")
    df <- utils::head(as.data.frame(df), max_rows)
    hdr <- paste0("<th>", html_escape(names(df)), "</th>", collapse = "")
    body <- apply(df, 1, function(r) {
      paste0("<tr>", paste0("<td>", html_escape(format(r, digits = 3)), "</td>", collapse = ""), "</tr>")
    })
    paste0("<table><thead><tr>", hdr, "</tr></thead><tbody>", paste(body, collapse = ""), "</tbody></table>")
  }

  top_de <- if (!is.null(de_table)) {
    d <- utils::head(de_table[de_table$significant, c("gene_name", "log2FoldChange", "padj", "direction")], 25)
    fmt_table(d)
  } else "<p><em>not run</em></p>"

  top_dmp <- if (!is.null(dmp_table)) {
    fmt_table(utils::head(dmp_table[dmp_table$significant, c("probe", "gene", "delta_beta", "padj", "direction")], 25))
  } else "<p><em>not run</em></p>"

  fig_html <- paste0(vapply(figs, function(f) {
    sprintf('<figure><img src="figures/%s" alt="%s"><figcaption>%s</figcaption></figure>', f, f, f)
  }, character(1)), collapse = "\n")

  tab_html <- paste0("<ul>", paste0(sprintf('<li><a href="tables/%s">%s</a></li>', tabs, tabs), collapse = ""), "</ul>")

  html <- sprintf('<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<style>
 :root{--fg:#1b1f23;--bg:#fff;--muted:#6a737d;--line:#e1e4e8;--accent:#0b5fff}
 @media (prefers-color-scheme:dark){:root{--fg:#e6edf3;--bg:#0d1117;--muted:#8b949e;--line:#30363d;--accent:#58a6ff}}
 body{margin:0 auto;max-width:1100px;padding:2rem 1.25rem;font:16px/1.6 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg)}
 h1{margin-bottom:.2rem} h2{margin-top:2.5rem;border-bottom:1px solid var(--line);padding-bottom:.3rem}
 .meta{color:var(--muted)} a{color:var(--accent)}
 table{border-collapse:collapse;width:100%%;font-size:14px;display:block;overflow-x:auto}
 th,td{border:1px solid var(--line);padding:.35rem .5rem;text-align:left;white-space:nowrap}
 th{background:rgba(127,127,127,.1)}
 figure{margin:0 0 2rem}
 img{max-width:100%%;height:auto;border:1px solid var(--line);border-radius:6px}
 figcaption{color:var(--muted);font-size:13px;margin-top:.35rem}
</style></head><body>
<h1>%s</h1>
<p class="meta">%s &middot; %s vs %s &middot; generated %s</p>
<h2>Cohort</h2>%s
<h2>Top differentially expressed genes</h2>%s
<h2>Top differentially methylated probes</h2>%s
<h2>Figures</h2>%s
<h2>Tables</h2>%s
<h2>Reproducibility</h2>
<p>Parameters: <a href="run_config.resolved.yml">run_config.resolved.yml</a> &middot;
Environment: <a href="sessionInfo.txt">sessionInfo.txt</a></p>
</body></html>',
    html_escape(cfg$contrast_name), html_escape(cfg$contrast_name),
    html_escape(cfg$project), html_escape(cfg$cohort$treatment), html_escape(cfg$cohort$reference),
    format(Sys.time(), "%Y-%m-%d %H:%M"),
    fmt_table(cohort_summary), top_de, top_dmp, fig_html, tab_html)

  path <- file.path(cfg$results_dir, "report.html")
  writeLines(html, path)
  log_ok("report written to ", path)
  invisible(path)
}
