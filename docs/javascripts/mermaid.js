// Initialize Mermaid with theme matching the current color scheme.
document.addEventListener("DOMContentLoaded", function () {
  var scheme = document.body.getAttribute("data-md-color-scheme");
  var theme = scheme === "slate" ? "dark" : "default";
  mermaid.initialize({ startOnLoad: true, theme: theme });
});
