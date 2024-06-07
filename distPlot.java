
import j4np.data.base.DataFrame;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.CompositeNode;
import j4np.hipo5.data.Event;
import j4np.hipo5.data.Node;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.io.HipoWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import java.util.List;
import java.io.File;
import java.io.IOException;
import java.lang.Math;

import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.H2F;
import twig.graphics.TGCanvas;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;
import java.io.PrintWriter;

import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;

/**
 *
 * @author tyson
 */
public class distPlot {

  public distPlot() {

  }

  public static int getSector(final double phi) {
    // shift in positive-phi direction by 3.5 sectors, result in [0,2*pi):
    final double phiShifted = Math.IEEEremainder(phi+Math.PI/6,2.*Math.PI)+Math.PI;
    // shifted sector number: 
    final int sectorShifted = (int)(phiShifted / (Math.PI/3)) + 1;
    // rotate back to proper sector:
    return sectorShifted<=3 ? sectorShifted-3+6 : sectorShifted-3;
  }

  public int[] hasMatchedAndInd(CompositeNode pt_test, double[] part, double[] matching_part, int charge, int sector){

    double Pz_o=part[2];
    int[] matched= new int[2];
    matched[0]=-1;
    matched[1]=-1;

    double bestRes=999;
    int bestInd=-1;

    for(int i = 0; i < pt_test.getRows(); i++){
      if(pt_test.getInt(2,i)==1){
        double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
        if(pt_test.getInt(1,i)==charge && resPz<0.2 && pt_test.getInt(2,i)==sector){
          if(resPz<bestRes){
            bestRes=resPz;
            bestInd=i;
          }
        }
      }
    }

    if(bestInd!=-1){
      matched[0]=pt_test.getInt(3,bestInd);
      matched[1]=bestInd;
      matching_part[0]=pt_test.getDouble(4,bestInd);
      matching_part[1]=pt_test.getDouble(5,bestInd);
      matching_part[2]=pt_test.getDouble(6,bestInd);
    } else{
      matching_part[0]=0;
      matching_part[1]=0;
      matching_part[2]=0;
    }

    return matched;
  }

  public int nPartsInSect(CompositeNode pt, int Sect){

    int n=0;
    for(int i = 0; i < pt.getRows(); i++){
      int s=pt.getInt(2,i);
      if(s==Sect){n++;}
    }
    return n;

  }

  public void printInsta(CompositeNode pt){
    System.out.println("Insta");
    for(int i=0;i<pt.getRows();i++){
      System.out.printf("sector = %2d, charge = %3d, pid = %3d, resp= %3.3f ", pt.getInt(2,i),pt.getInt(1,i),pt.getInt(3,i),pt.getDouble(10,i));
      System.out.printf("p %9.5f  %9.5f %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i));
      System.out.printf("pcal ADC = %9.5f, ecin ADC = %9.5f, ecout ADC = %9.5f, sum HTCC ADC %9.5f \n", pt.getDouble(12,i),pt.getDouble(13,i),pt.getDouble(14,i),pt.getDouble(15,i));
      System.out.printf("HTCC PMTS Left: ");
      for(int k = 16; k < 20; k++) System.out.printf(" %9.5f ",pt.getDouble(k,i));
      System.out.printf("\nHTCC PMTS Right: ");
      for(int k = 20; k < 24; k++) System.out.printf(" %9.5f ",pt.getDouble(k,i));
      System.out.println("\n");
    }
  }
  public void printREC(CompositeNode pt_test, CompositeNode cf_test){
    System.out.println("REC");
    for(int i=0;i<pt_test.getRows();i++){
      System.out.printf("pindex = %2d, sector = %2d, charge = %3d, pid = %3d, status = %3d, track chi^2 = %9.5f ",cf_test.getInt(0, i), pt_test.getInt(2,i),pt_test.getInt(1,i),pt_test.getInt(3,i),pt_test.getInt(11,i),pt_test.getDouble(12,i));
      System.out.printf("p %9.5f  %9.5f %9.5f\n",pt_test.getDouble(4,i),pt_test.getDouble(5,i),pt_test.getDouble(6,i));
      System.out.printf("pcal e = %9.5f, ecin e = %9.5f, ecout e = %9.5f, nb nphe = %9.5f \n\n", pt_test.getDouble(13,i), pt_test.getDouble(14,i), pt_test.getDouble(15,i), pt_test.getDouble(16,i));
    }
  }

  public void printTraj(Bank rectraj,int det, double htcctheta,double htccphi){
    System.out.println("REC::Traj (at HTCC)");
    for(int i=0;i<rectraj.getRows();i++){
      if(rectraj.getInt("detector",i)==det){
        int pindex = rectraj.getInt("pindex",i);
        double x = rectraj.getFloat("x",i);
        double y = rectraj.getFloat("y",i);
        double z = rectraj.getFloat("z",i);
        double theta = Math.atan2(Math.sqrt(x*x+y*y),z);
        double phi   = Math.atan2(y,x);

        double dtheta=(htcctheta-theta)*180/Math.PI;
        //Wrap delta-phi into -pi/pi as in REC
        double dphi=Math.IEEEremainder(htccphi-phi,2.*Math.PI)*180/Math.PI;
        System.out.printf("pindex %d x %.3f y %.3f z %.3f theta %.3f phi %.3f dtheta %.3f dphi %.3f \n",pindex,x,y,z,theta*180/Math.PI,phi*180/Math.PI,dtheta,dphi);
      }
    }

    System.out.println("\n");
    //rectraj.show();
  }

  public int findTrajRow(Bank rectraj,int det, int pindex){
    int row =-1;
    for(int i=0;i<rectraj.getRows();i++){
      if(rectraj.getInt("detector",i)==det && rectraj.getInt("pindex",i)==pindex){
        row=i;
      }
    }
    return row;
  }

  public void plot(String file, int nEv){

    HipoReader r = new HipoReader(file); //_noRecEl_InstaEl_
    System.out.println("read "+file);
    CompositeNode tr = new CompositeNode(32100,2,"i",1200);
    CompositeNode pt = new CompositeNode(32100,3,"i",1200);
    CompositeNode cf = new CompositeNode(32100,11,"i",1200);
    CompositeNode pt_test = new CompositeNode(32100,22,"i",1200);
    CompositeNode cf_test = new CompositeNode(32100,21,"i",1200);
    Bank htccrec = r.getBank("HTCC::rec");
    Bank htccadc = r.getBank("HTCC::adc");
    Bank rectraj = r.getBank("RECAI::Traj");
    Event event = new Event();
    int counter = 0;

    H2F hThetaPhi = new H2F("dPhi vs dTheta",40,-20,20,72,-36,36);
    hThetaPhi.attr().setTitleX("dTheta [Degrees]");
    hThetaPhi.attr().setTitleY("dPhi [Degrees]");

    H2F hThetaPhi_notEl = new H2F("dPhi vs dTheta",40,-20,20,72,-36,36);
    hThetaPhi_notEl.attr().setTitleX("dTheta [Degrees]");
    hThetaPhi_notEl.attr().setTitleY("dPhi [Degrees]");

    H2F hThetaPhi_m = new H2F("dPhi vs dTheta",40,-20,20,72,-36,36);
    hThetaPhi_m.attr().setTitleX("dTheta [Degrees]");
    hThetaPhi_m.attr().setTitleY("dPhi [Degrees]");

    while(r.hasNext() && counter<nEv){
      r.next(event);
      event.read(tr,32100,2);
      event.read(pt,32100,3);
      event.read(cf,32100,11);
      event.read(pt_test,32100,22);
      event.read(cf_test,32100,21);
      event.read(htccrec);
      event.read(htccadc);
      event.read(rectraj);
      int nrows = tr.getRows();
      
      int shown=0,isMatched=-1;
      
      int nrowstest=pt_test.getRows();
      //for simplicity only want one particle in sector
      //use sector 1 as already skimmed to check if have electron in sector 1
      if(nPartsInSect(pt, 1)==1 && nPartsInSect(pt_test, 1)==1){
        for(int i = 0; i < nrowstest; i++){
          if(pt_test.getInt(2,i)==1 && pt_test.getInt(1,i)==-1){
            double Px_o=pt_test.getDouble(4,i);
            double Py_o=pt_test.getDouble(5,i);
            double Pz_o=pt_test.getDouble(6,i);
            double[] part= new double[3];
            part[0]=Px_o;
            part[1]=Py_o;
            part[2]=Pz_o;
            double[] match= new double[3];
            int[] isMatchedArr=hasMatchedAndInd(pt, part, match, -1,1);
            //have a reconstructed cluster in HTCC
            //only take one cluster in HTCC for simplicity
            int trajrow = findTrajRow(rectraj, 15, cf_test.getInt(0, i));
            if(htccrec.getRows()==1 && trajrow!=-1){

              double x = rectraj.getFloat("x",trajrow);
              double y = rectraj.getFloat("y",trajrow);
              double z = rectraj.getFloat("z",trajrow);
              double theta = Math.atan2(Math.sqrt(x*x+y*y),z);
              double phi   = Math.atan2(y,x);

              double htcctheta = htccrec.getFloat("theta", 0);
              double htccphi = htccrec.getFloat("phi", 0);
              double htcc_sector=getSector(htccphi);
              /*double htccx = htccrec.getFloat("x",0);
              double htccy = htccrec.getFloat("y",0);
              double htccz = htccrec.getFloat("z",0);
              double htcctheta = Math.atan2(Math.sqrt(htccx*htccx+htccy*htccy),htccz);
              double htccphi   = Math.atan2(htccy,htccx);*/


              double dtheta=(htcctheta-theta)*180/Math.PI;
              //Wrap delta-phi into -pi/pi as in REC
              double dphi=Math.IEEEremainder(htccphi-phi,2.*Math.PI)*180/Math.PI;
              
              //part is matched to InstaRec electron
              if(isMatchedArr[0]==11 ){
                
                //part isn't IDed as electron as HTCC cluster not assigned to part
                if(pt_test.getInt(3,i)!=11 && pt_test.getDouble(16,i)==0){
                  hThetaPhi_m.fill(dtheta,dphi);

                  if(Math.abs(dtheta)<10 && Math.abs(dphi)<18){
                    System.out.printf("Insta 11 REC !11 but dTheta %.3f & dPhi %.3f \n",dtheta,dphi);
                    printInsta(pt);
                    printREC(pt_test,cf_test);
                    System.out.printf("HTC::rec theta %f phi %f sector %f \n",htcctheta*180/Math.PI,htccphi*180/Math.PI,htcc_sector);
                    htccrec.show();
                    System.out.println("HTC::adc");
                    htccadc.show();
                    printTraj(rectraj,15,htcctheta,htccphi);
                  }

                }
                //part is an electron, must have HTCC cluster assigned
                else {
                  hThetaPhi.fill(dtheta,dphi);

                  if(Math.abs(dtheta)>=10 && Math.abs(dphi)>=18){
                    System.out.printf("Insta 11 REC 11 but dTheta %.3f & dPhi %.3f \n",dtheta,dphi);
                    printInsta(pt);
                    printREC(pt_test,cf_test);
                  }

                }
                
              }
              //part is matched to InstaRec electron
              else{
                
                //part isn't IDed as electron as HTCC cluster not assigned to part
                if(pt_test.getInt(3,i)!=11 && pt_test.getDouble(16,i)==0){
                  hThetaPhi_notEl.fill(dtheta,dphi);
                }
                
              }

            }

          }
        }

      }
      

    }

    TGCanvas cd = new TGCanvas();
    cd.setTitle("InstaRec 11, REC !11");
    cd.draw(hThetaPhi_m);

    TGCanvas cd3 = new TGCanvas();
    cd3.setTitle("InstaRec !11, REC !11");
    cd3.draw(hThetaPhi_notEl);

    TGCanvas cd2 = new TGCanvas();
    cd2.setTitle("InstaRec 11, REC 11");
    cd2.draw(hThetaPhi);
      

  }
  
  public static void main(String[] args) {

    //String field="inbending";
    String field="outbending";
    int nEvs=100000;

    String file="output_noRecEl_InstaEl_"+field+"_RfRecooked.h5";

    distPlot pter = new distPlot();
    pter.plot(file,nEvs);

  }

}