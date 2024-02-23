import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.base.Detector;
import org.jlab.geom.component.ScintillatorPaddle;
import org.jlab.geom.prim.Point3D;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;
import java.io.PrintWriter;

//to run ~/coatjava/bin/jshell_r.sh
// /open gen_L_strip_convTable.java
// genTable()

void genTable(){
    Detector d = GeometryFactory.getDetector(DetectorType.ECAL);

    //d.getSector(N) -- gives sector N
    //distance between strips doesn't change by sector
    //d.getSector(0).getSuperlayer(N) -- gives PCAL at 0, ECIN at 1, ECOUT at 3
    //d.getSector(0).getSuperlayer(0).getLayer(N) -- gives N layers
    //layer 0 is U view on front facing block, 1 is V View, 2 is W View
    //there's N blocks of 3 Views, only care about the front facing one as this is where distances are reported
    //d.getSector(0).getSuperlayer(0).getLayer(0).getComponent(N) -- gives the Nth strip

    Path filePath = Paths.get("LtoStrip_convTable.csv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath))) {
        for(int det=0;det<3;det++){
            for (int view=0;view<3;view++){

                StringBuilder csvLineBuilder = new StringBuilder();
                
                //wasn't sure where the length is calculated from
                //it's calculated from front facing strip ie start=0
                int start=0;
                /*int start=6;
                if(det==2){
                    start=12;
                }*/
                /*int start=12;
                if(det==2){
                    start=21;
                }*/

                Point3D p0 = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(0)).getLine().origin();
                if(det==0 && view>0){
                    //from reco in ECCommon
                    // The distance from the edge on PCAL
                    // for layers V and W are calculated from the wider
                    // strips edge.
                    p0 = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(61)).getLine().origin();
                } else{
                    // for PCAL V, W don't want to start with 0
                    csvLineBuilder.append(String.format("%.6f,", 0.0));
                }


                
                //not all layers have same number of components
                //to make reading table easier, just pad the later ones with 0s
                //int nComp=d.getSector(0).getSuperlayer(det).getLayer(view).getNumComponents();
                for(int comp=1;comp<69;comp++){
                    double dist=0.0;
                    if(det==0){
                        //68 strips in PCAL U view
                        if(view==0){
                            if(comp<68){
                                Point3D pN = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp)).getLine().origin();
                                dist=p0.distance(pN); 
                            } else{
                                //distance to end of last strip?
                                //take distance to start of strip + length of strip before
                                Point3D pN = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp-1)).getLine().origin();
                                Point3D pNe = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp-2)).getLine().origin();
                                dist=p0.distance(pN)+pN.distance(pNe);
                            }
                            
                        } else{
                            //62 strips in PCAL V, W views
                            if (comp==0){
                                //distance to end of last strip?
                                //take distance to start of strip + length of strip before
                                //last strip comes first for pcal V W
                                Point3D pN = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp-1)).getLine().origin();
                                Point3D pNe = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp-2)).getLine().origin();
                                double hL = 394.2*0.5;
                                double hyp = Math.sqrt(hL*hL + 385.2*385.2);
                                double theta = Math.acos(hL/hyp);
                                double proj  = 4.5*Math.cos(theta);
                                dist=p0.distance(pN)+pN.distance(pNe)+proj;
                            } else if(comp<62){
                                //from reco in ECCommon
                                // The distance from the edge on PCAL
                                // for layers V and W are calculated from the wider
                                // strips edge.
                                Point3D pN = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp)).getLine().origin();
                                double hL = 394.2*0.5;
                                double hyp = Math.sqrt(hL*hL + 385.2*385.2);
                                double theta = Math.acos(hL/hyp);
                                double proj  = 4.5*Math.cos(theta);
                                dist=p0.distance(pN)+proj;

                            } 
                        }
                    } else{
                        //36 strips U/V/W ECIN, ECOUT
                        if(comp<36){
                            Point3D pN = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp)).getLine().origin();
                            dist=p0.distance(pN);
                        } else if (comp==36){
                            //distance to end of last strip?
                            //take distance to start of strip + length of strip before
                            Point3D pN = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp-1)).getLine().origin();
                            Point3D pNe = ((ScintillatorPaddle) d.getSector(0).getSuperlayer(det).getLayer(start+view).getComponent(comp-2)).getLine().origin();
                            dist=p0.distance(pN)+pN.distance(pNe);
                        }
                    }
                    csvLineBuilder.append(String.format("%.6f,", dist));
                }

                if(det==0 && view>0){
                    // for PCAL V, W add extra at end 0
                    csvLineBuilder.append(String.format("%.6f,", 0.0));
                }

                csvLineBuilder.deleteCharAt(csvLineBuilder.length() - 1);
                writer.println(csvLineBuilder.toString());
            }
        }
    }catch (IOException e) {
        e.printStackTrace();
    }

}


