pandoc -s --filter pandoc-crossref -M "crossrefYaml=crossref_config.yaml" --template=github.html --mathjax --highlight-style=tango README.md -o index.html
