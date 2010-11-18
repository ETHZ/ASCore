void rootlogon() {
  gSystem->SetIncludePath(" -Iinclude/");
  gSystem->Load("libPhysics");
  gSystem->Load("shlib/libDiLeptonAnalysis.so");
}
