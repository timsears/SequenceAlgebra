with import <nixpkgs> {};
let hp = haskell.packages.ghc844;
in
stdenv.mkDerivation {
  name = "SequenceAlgebra";
  buildInputs = [
  (hp.ghcWithHoogle (p : [
                          #none yet
                           ] ) )
    hp.cabal-install
  ];
  src = null;
}
