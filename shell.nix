with import <nixpkgs> {};
let hp = haskell.packages.ghc844;
in
stdenv.mkDerivation {
  name = "SequenceAlgebra";
  buildInputs = [
  (hp.ghcWithHoogle (p : [p.QuickCheck p.optparse-generic p.boltzmann-samplers] ) )
    hp.cabal-install
  ];
  src = null;
}
