#!/bin/sh
# Assembles variants.html from the _p*.html parts and the two SVG figures.
#
# The XML prolog is stripped so the SVG can be inlined, and figure 2 gets its
# glyph ids renamed: both figures carry <use xlink:href="#glyph-N-M"> and in one
# document the second would otherwise resolve every glyph against the first.
set -e
cd "$(dirname "$0")"

tail -n +2 fig-coverage.svg                        > _figA.svg
tail -n +2 fig-lipsey-vs-n.svg | sed 's/glyph-/gB-/g' > _figB.svg

cat _p0.html _p1.html _p2.html _p3.html \
    _p4.html _figA.svg \
    _p5.html _figA.svg \
    _p6.html _figA.svg _p7.html _figB.svg \
    _p8.html > variants.html

echo "wrote variants.html ($(wc -c < variants.html) bytes)"
