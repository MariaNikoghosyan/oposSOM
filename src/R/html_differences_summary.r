pipeline.htmlDifferencesSummary <- function()
{
  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Differences",
                        "0verview.html")

  util.info("Writing:", filename)
  f <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Differences Summary of ", files.name, " dataset</title>
    <style>
      body {
        margin: 0;
        padding: 0;
        color: #333;
        background-color: #fff;
        font: normal normal normal 14px/180% sans-serif;
      }

      #wrapper {
        width: 90%;
        min-width: 400px;
        max-width: 800px;
        margin: 20px auto;
      }

      h1, h2 {
        margin: 30px 0 0 0;
        line-height: 210%;
        border-bottom: 1px solid #eee;
      }

      dl {
        line-height: 180%;
      }

      dl dt {
        width: 50%;
        float: left;
        color: #111;
      }

      dl dt:after {
        content: ':';
      }

      a {
        color: #4183C4;
        text-decoration: none;
      }

      table {
        width: 100%;
        margin: 24px 0;
        border-collapse: collapse;
      }

      table th,
      table td {
        text-align: left;
        padding: 4px 6px;
      }

      table thead tr,
      table tbody tr:nth-child(2n) {
        background-color: #f0f0f0;
      }

      table a {
        display: block;
      }
    </style>
  </head>
  <body>
    <div id=\"wrapper\">
      <h1>Differences Maps Summary Sheets</h1>

      <p>
        Some fancy text here...
      </p>

      <h2>Something</h2>

      <p>
        Another fancy text here...
      </p>

      <table>
        <thead>
          <tr>
            <th>Categories</th>
            <th>Summary Sheets</th>
            <th>Global Gene List</th>
            <th>Local Gene Lists</th>
            <th>Gene Set List</th>
          </tr>
        </thead>
        <tbody>", sep="", file=f)

  for (m in 1:ncol(indata))
  {
    name <- colnames(indata)[m]
    fname <- make.names(name)

    cat("
          <tr>
            <td>", name, "</td>
            <td>
              <a href=\"Reports/", fname, ".pdf\" target=\"_blank\">
                PDF
              </a>
            </td>
            <td>
              <a href=\"CSV Sheets/Gene Lists - Global/", fname, ".csv\" target=\"_blank\">
                CSV
              </a>
            </td>
            <td>", sep="", file=f)

    for (i in 1:length(GS.infos.samples[[m]]$spots))
    {
      cat("
              <a href=\"CSV Sheets/Gene Lists - Local/", fname, ".", i, ".csv\" target=\"_blank\">
                CSV ", i, "
              </a>", sep="", file=f)
    }

    cat("
            </td>
            <td>
              <a href=\"CSV Sheets/Gene Set Lists - Global/", fname, ".csv\" target=\"_blank\">
                CSV
              </a>
            </td>
          </tr>", sep="", file=f)
  }

  cat("
        </tbody>
      </table>
    </div>
  </body>
</html>", sep="", file=f)

  close(f)
}
