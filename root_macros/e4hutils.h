auto getPhi(double x, double y, double z){return ROOT::Math::XYZVector{x, y, z}.Phi();}
auto getTheta(double x, double y, double z){return ROOT::Math::XYZVector{x, y, z}.Theta();}
auto getEnergy(double px, double py, double pz, double mass){return ROOT::Math::PxPyPzMVector{px, py, pz, mass}.E();}
auto getPt(double px, double py, double pz, double mass){return ROOT::Math::PxPyPzMVector{px, py, pz, mass}.Pt();}
