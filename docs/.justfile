list:
  just --list

build:
  julia --color=yes --project make.jl

serve:
  julia --project -e 'using LiveServer; serve(dir="html")'

clean:
  rm -rf build html pdf
  rm -f src/figs/*.html src/figs/*.pdf
