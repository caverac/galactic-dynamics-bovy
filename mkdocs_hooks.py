"""MkDocs build hooks.

The mkdocs-bibtex plugin renders citations through an APA CSL style whose
output contains pandoc-style no-case spans, e.g. ``[Xue et al.]{.nocase}``.
Python-Markdown does not understand that pandoc attribute syntax, so the markup
leaks into the page as literal text. This hook strips the span wrappers,
leaving just the enclosed text, for every rendered page.
"""

import re
from typing import Any

_NOCASE_SPAN = re.compile(r"\[([^\[\]]+)\]\{\.nocase\}")


def on_page_content(html: str, **_kwargs: Any) -> str:
    """Remove leftover pandoc no-case spans from rendered HTML.

    Parameters
    ----------
    html : str
        The page HTML produced from Markdown.
    **_kwargs
        Other MkDocs event arguments (page, config, files); unused.

    Returns
    -------
    str
        The HTML with ``[text]{.nocase}`` wrappers reduced to ``text``.
    """
    return _NOCASE_SPAN.sub(r"\1", html)
