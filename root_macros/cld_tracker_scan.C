int cld_tracker_scan(){
  //events->Draw("sqrt(OuterTrackerEndcapCollection.position.x*OuterTrackerEndcapCollection.position.x + OuterTrackerEndcapCollection.position.y*OuterTrackerEndcapCollection.position.y):abs(OuterTrackerEndcapCollection.position.z)");
  
  //events->Draw("sqrt(OuterTrackerEndcapCollection.position.x*OuterTrackerEndcapCollection.position.x + OuterTrackerEndcapCollection.position.y*OuterTrackerEndcapCollection.position.y):OuterTrackerEndcapCollection.position.z");
  //events->Draw("sqrt(OuterTrackerBarrelCollection.position.x*OuterTrackerBarrelCollection.position.x + OuterTrackerBarrelCollection.position.y*OuterTrackerBarrelCollection.position.y):OuterTrackerBarrelCollection.position.z", "", "same");

  events->Draw("OuterTrackerEndcapCollection.position.x:OuterTrackerEndcapCollection.position.z");
  events->Draw("sqrt(OuterTrackerBarrelCollection.position.x*OuterTrackerBarrelCollection.position.x + OuterTrackerBarrelCollection.position.y*OuterTrackerBarrelCollection.position.y):OuterTrackerBarrelCollection.position.z", "", "same");

  //events->Draw("InnerTrackerEndcapCollection.position.x:InnerTrackerEndcapCollection.position.z", "", "same");
  events->Draw("sqrt(InnerTrackerBarrelCollection.position.x*InnerTrackerBarrelCollection.position.x + InnerTrackerBarrelCollection.position.y*InnerTrackerBarrelCollection.position.y):InnerTrackerBarrelCollection.position.z", "", "same");

  events->Draw("VertexEndcapCollection.position.x:VertexEndcapCollection.position.z", "", "same");
  events->Draw("sqrt(VertexBarrelCollection.position.x*VertexBarrelCollection.position.x + VertexBarrelCollection.position.y*VertexBarrelCollection.position.y):VertexBarrelCollection.position.z", "", "same");

  TLine* l = new TLine(0,0,3000*cos(0.872665),3000*sin(0.872665));
  l->Draw("same");
  return true;
}
