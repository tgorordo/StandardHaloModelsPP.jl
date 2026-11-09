{
  description = "StandardHaloModelsPP.jl development and documentation environment";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
    in
    {
      devShells = forAllSystems (
        pkgs:
        let
          # Dependencies for docs/make.jl esp. the LaTeX/PDF build and the Asymptote figs.
          # Documenter's documenter.sty loads fontspec/polyglossia (luatex), minted,
          # tabulary, xltabular, tcolorbox, adjustbox, booktabs, and sets DejaVu fonts.
          texlive = pkgs.texliveMedium.withPackages (
            ps: with ps; [
              collection-latexextra
              collection-luatex
              collection-fontsrecommended
              latexmk
              dejavu
              starfont
            ]
          );
        in
        {
          default = pkgs.mkShell {
            packages = [
              pkgs.julia-bin
              pkgs.just

              texlive # also provides asy (Asymptote) and latexminted
              pkgs.ghostscript
              (pkgs.python3.withPackages (ps: [ ps.pygments ])) # minted highlighting
            ];

            # DejaVu via fontconfig, for anything resolving fonts outside kpathsea.
            FONTCONFIG_FILE = pkgs.makeFontsConf { fontDirectories = [ pkgs.dejavu_fonts ]; };
          };
        }
      );

      formatter = forAllSystems (pkgs: pkgs.nixfmt);
    };
}
