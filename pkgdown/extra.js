/* The articles as a book.

   The chapters, their order and the parts are read from the navbar, which
   pkgdown builds from `articles:` in _pkgdown.yml: chapter 1 is "Get started",
   then every article of the Articles menu in order, grouped by the menu's
   headings (the parts). Adding an article to `articles:` adds a chapter and
   renumbers the ones after it; nothing here needs to change.

   On a chapter page this script
   - numbers the title and the section headings (5, 5.1, 5.1.1),
   - shows the chapters in the right sidebar in place of "On this page", with
     the current chapter's sections nested under it and the section in view
     highlighted (on a phone, a collapsed "Chapters" box below the title),
   - adds links to the previous and the next chapter at the bottom.
   On the articles index it numbers the entries the same way. */
(function () {
  "use strict";

  function pathOf(href) {
    var p = new URL(href, document.baseURI).pathname;
    return p.replace(/\/index\.html$/, "/").replace(/\.html$/, "");
  }

  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) {
      if (k === "text") e.textContent = attrs[k];
      else e.setAttribute(k, attrs[k]);
    });
    (children || []).forEach(function (c) { if (c) e.appendChild(c); });
    return e;
  }

  /* The book: [{number, title, href, path, part}] in reading order. */
  function readBook() {
    var chapters = [];
    var add = function (a, part) {
      if (/\/articles\/(index\.html)?$/.test(new URL(a.href).pathname)) return;
      chapters.push({ title: a.textContent.trim(), href: a.href, path: pathOf(a.href), part: part });
    };
    document.querySelectorAll("#navbar .navbar-nav > li > a.nav-link").forEach(function (a) {
      if (/\/articles\/[^/]+\.html$/.test(new URL(a.href).pathname)) add(a, null);
    });
    var menu = document.querySelector('#dropdown-articles + .dropdown-menu');
    var part = null;
    if (menu) {
      menu.querySelectorAll(".dropdown-header, a.dropdown-item").forEach(function (x) {
        if (x.classList.contains("dropdown-header")) part = x.textContent.trim();
        else add(x, part);
      });
    }
    chapters.forEach(function (c, i) { c.number = i + 1; });
    return chapters;
  }

  /* Number the section headings of a chapter page (or of a fetched one) and
     return them as a tree [{number, title, id, children}]. */
  function numberSections(root, chapter, edit) {
    var tree = [];
    var walk = function (container, level, prefix, into) {
      var k = 0;
      container.querySelectorAll(":scope > .section.level" + level).forEach(function (sec) {
        var h = sec.querySelector(":scope > h" + level);
        if (!h || h.hasAttribute("data-toc-skip")) return;
        k += 1;
        var number = prefix + "." + k;
        var copy = h.cloneNode(true);
        copy.querySelectorAll(".anchor, .section-number").forEach(function (x) { x.remove(); });
        var node = { number: number, title: copy.innerHTML.trim(), id: sec.id || h.id, heading: h, children: [] };
        if (edit) h.insertBefore(el("span", { class: "section-number", text: number }), h.firstChild);
        into.push(node);
        if (level < 3) walk(sec, level + 1, number, node.children);
      });
    };
    var main = root.querySelector("main");
    if (main) {
      // the sections are children of main, or of a wrapper div
      var top = main.querySelector(".section.level2");
      if (top) walk(top.parentElement, 2, String(chapter.number), tree);
    }
    return tree;
  }

  function sectionList(nodes, href) {
    if (!nodes.length) return null;
    return el("ol", { class: "book-sections" }, nodes.map(function (n) {
      var a = el("a", { href: href + "#" + n.id, "data-target": n.id }, [
        el("span", { class: "book-num", text: n.number }),
        el("span", { class: "book-title" })
      ]);
      a.lastChild.innerHTML = n.title; // the heading's own markup, e.g. <code>
      return el("li", {}, [a, sectionList(n.children, href)]);
    }));
  }

  function chapterLink(c, current) {
    var a = el("a", { href: c.href, class: "book-chapter-link" }, [
      el("span", { class: "book-num", text: String(c.number) }),
      el("span", { class: "book-title", text: c.title })
    ]);
    if (current) a.setAttribute("aria-current", "page");
    return a;
  }

  /* The other chapters' sections are read from their pages when asked for. */
  var fetched = {};
  function unfoldButton(c, li) {
    var b = el("button", { type: "button", class: "book-unfold", "aria-expanded": "false",
                           "aria-label": "Sections of chapter " + c.number });
    b.addEventListener("click", function () {
      var open = b.getAttribute("aria-expanded") === "true";
      var list = li.querySelector(":scope > .book-sections");
      if (open) {
        b.setAttribute("aria-expanded", "false");
        if (list) list.hidden = true;
        return;
      }
      b.setAttribute("aria-expanded", "true");
      if (list) { list.hidden = false; return; }
      (fetched[c.path] = fetched[c.path] || fetch(c.href).then(function (r) {
        if (!r.ok) throw new Error(r.status);
        return r.text();
      }).then(function (html) {
        var doc = new DOMParser().parseFromString(html, "text/html");
        return numberSections(doc, c, false);
      })).then(function (tree) {
        var ol = sectionList(tree, c.href);
        if (ol) li.appendChild(ol);
      }).catch(function () {
        b.hidden = true;
      });
    });
    return b;
  }

  function buildPanel(book, current, tree, idSuffix) {
    var nav = el("nav", { class: "book-panel", "aria-label": "Chapters" });
    var list = el("ol", { class: "book-chapters" });
    var canFetch = location.protocol !== "file:";
    var part;
    book.forEach(function (c) {
      if (c.part !== null && c.part !== part) {
        list.appendChild(el("li", { class: "book-part", text: c.part }));
      }
      part = c.part;
      var isCurrent = current && c.path === current.path;
      var li = el("li", { class: "book-chapter" + (isCurrent ? " current" : "") });
      var row = el("div", { class: "book-row" }, [chapterLink(c, isCurrent)]);
      li.appendChild(row);
      if (isCurrent) {
        var ol = sectionList(tree, "");
        if (ol) li.appendChild(ol);
      } else if (canFetch) {
        row.appendChild(unfoldButton(c, li));
      }
      list.appendChild(li);
    });
    nav.appendChild(el("h2", { text: "Chapters" }));
    nav.appendChild(list);
    return nav;
  }

  function pager(book, i) {
    var prev = book[i - 1], next = book[i + 1];
    var side = function (c, rel, label) {
      if (!c) return el("span");
      return el("a", { href: c.href, rel: rel, class: "book-pager-" + rel }, [
        el("span", { class: "book-pager-label", text: label }),
        el("span", { class: "book-pager-title", text: c.number + " " + c.title })
      ]);
    };
    return el("nav", { class: "book-pager", "aria-label": "Previous and next chapter" }, [
      side(prev, "prev", "← Previous chapter"),
      side(next, "next", "Next chapter →")
    ]);
  }

  /* Highlight the section in view, in every panel on the page. */
  function spy(tree, panels) {
    var flat = [];
    (function f(ns) { ns.forEach(function (n) { flat.push(n); f(n.children); }); })(tree);
    if (!flat.length) return;
    var last = null;
    var update = function () {
      var cut = Math.max(100, window.innerHeight * 0.25);
      var active = null;
      flat.forEach(function (n) {
        if (n.heading.getBoundingClientRect().top <= cut) active = n;
      });
      // at the end of the page, the last section counts as in view
      var end = document.documentElement;
      if (window.innerHeight + window.scrollY >= end.scrollHeight - 2) active = flat[flat.length - 1];
      var id = active ? active.id : null;
      if (id === last) return;
      last = id;
      panels.forEach(function (p) {
        p.querySelectorAll("a[data-target]").forEach(function (a) {
          var on = a.getAttribute("data-target") === id;
          a.classList.toggle("active", on);
          if (on && p.offsetParent !== null && p.scrollHeight > p.clientHeight) {
            var top = a.offsetTop - p.offsetTop;
            if (top < p.scrollTop || top > p.scrollTop + p.clientHeight - 40) {
              p.scrollTop = top - p.clientHeight / 3;
            }
          }
        });
      });
    };
    var queued = false;
    window.addEventListener("scroll", function () {
      if (queued) return;
      queued = true;
      requestAnimationFrame(function () { queued = false; update(); });
    }, { passive: true });
    update();
  }

  function chapterPage(book, i) {
    var c = book[i];
    var main = document.querySelector("main#main");
    var header = main.querySelector(".page-header");
    var h1 = header && header.querySelector("h1");
    if (h1) h1.insertBefore(el("span", { class: "section-number", text: String(c.number) }), h1.firstChild);
    var tree = numberSections(document, c, true);

    var aside = main.parentElement.querySelector(":scope > aside");
    var panels = [];
    if (aside) {
      var toc = aside.querySelector("#toc");
      if (toc) toc.hidden = true;
      var side = buildPanel(book, c, tree);
      side.classList.add("book-panel-side");
      aside.insertBefore(side, aside.firstChild);
      aside.classList.add("book-aside");
      panels.push(side);
    }
    var box = el("details", { class: "book-panel-phone" }, [el("summary", { text: "Chapters" })]);
    var phone = buildPanel(book, c, tree);
    box.appendChild(phone);
    box.addEventListener("click", function (e) {
      if (e.target.closest && e.target.closest("a")) box.open = false;
    });
    if (header) header.after(box);
    panels.push(phone);

    main.appendChild(pager(book, i));
    spy(tree, panels);
  }

  function indexPage(book) {
    var byPath = {};
    book.forEach(function (c) { byPath[c.path] = c; });
    document.querySelectorAll("main#main dl dt a").forEach(function (a) {
      var c = byPath[pathOf(a.href)];
      if (!c) return;
      var dt = a.closest("dt");
      dt.insertBefore(el("span", { class: "section-number", text: String(c.number) }), dt.firstChild);
      dt.classList.add("book-entry");
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    var book = readBook();
    if (!book.length) return;
    var here = pathOf(location.href);
    for (var i = 0; i < book.length; i++) {
      if (book[i].path === here) return chapterPage(book, i);
    }
    if (document.querySelector(".template-article-index")) indexPage(book);
  });
})();
