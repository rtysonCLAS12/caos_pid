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
import twig.data.StatNumber;
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
public class plotter {

  public plotter() {

  }

  public int hasMatched(CompositeNode pt_test, double[] part, double[] matching_part, int charge){

    double Pz_o=part[2];
    int matched=-1;

    double bestRes=999;
    int bestInd=-1;

    for(int i = 0; i < pt_test.getRows(); i++){
      if(pt_test.getInt(2,i)==1){
        double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
        if(pt_test.getInt(1,i)==charge && resPz<0.2){
          if(resPz<bestRes){
            bestRes=resPz;
            bestInd=i;
          }
        }
      }
    }

    if(bestInd!=-1){
      matched=pt_test.getInt(3,bestInd);
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


  public void plot(String file,String varName, String varUnits, double low, double high, double step, Boolean requireMatching){

    GraphErrors gEff = new GraphErrors();
    gEff.attr().setMarkerColor(2);
    gEff.attr().setMarkerSize(10);
    gEff.attr().setTitle("Level3 Efficiency");
    gEff.attr().setTitleX(varName + " " + varUnits);
    gEff.attr().setTitleY("Metrics");

    GraphErrors gPur = new GraphErrors();
    gPur.attr().setMarkerColor(5);
    gPur.attr().setMarkerSize(10);
    gPur.attr().setTitle("Level3 Purity");
    gPur.attr().setTitleX(varName + " " + varUnits);
    gPur.attr().setTitleY("Metrics");


    for (double var = low; var < high; var += step) {
      //doing weird thing so i can mimick HTCC binning
      if(varName=="Theta" && var==21){step=7.0;}

      HipoReader r = new HipoReader(file);

      CompositeNode pt = new CompositeNode(32100,3,"sssifffffffsffffffffffff",1200);
      CompositeNode pt_test = new CompositeNode(32100,22,"sssifffffffsfffff",1200);

      Event event = new Event();
      int counter = 0;

      int hasInstaEl_hasAIEl=0, hasInstaEl=0;
      int hasAIEl=0, hasAIEl_hasInstaEl=0;


      while(r.hasNext()){
          counter++;
          r.next(event);

          //System.out.println("reading");

          event.read(pt,32100,3);
          event.read(pt_test,32100,22);
          int nrows = pt.getRows();
          int nrowstest=pt_test.getRows();
          //int nrowstest_nAI=pt_test_nAI.getRows();

          //System.out.printf(" nrows insta %d nrows ai rec %d \n",nrows,nrowstest);

          //check if InstaRec has e-
          for(int j = 0; j < nrows; j++){

            if(pt.getInt(2,j)!=0){
              //has Insta El
              if(pt.getInt(3,j) == 11 ){//&& pt.getDouble(15,j)>0.003

                double Px_o=pt.getDouble(4,j);
                double Py_o=pt.getDouble(5,j);
                double Pz_o=pt.getDouble(6,j);
                double P_o=Math.sqrt(Px_o*Px_o+Py_o*Py_o+Pz_o*Pz_o);
                double Theta_o = (float) Math.acos((float)Pz_o / P_o)*(180.0 / Math.PI);// Math.atan2(Math.sqrt(Px*Px+Py*Py),Pz);
                double Phi_o = (float) Math.atan2(Py_o, Px_o)*(180.0 / Math.PI);

                int isMatched=1;
                if(requireMatching){
                  double[] el= new double[3];
                  el[0]=Px_o;
                  el[1]=Py_o;
                  el[2]=Pz_o;
                  double[] match= new double[3];
                  isMatched=hasMatched(pt_test, el, match, -1);
                }

                double var_to_check=P_o;
                if(varName=="Theta"){
                  var_to_check=Theta_o;
                } else if(varName=="Phi"){
                  var_to_check=Phi_o;
                }

                if(isMatched!=-1){
                  if(var_to_check>=var && var_to_check<(var+step)){
                    hasInstaEl++;
  
                    for(int i = 0; i < nrowstest; i++){
                      if(pt_test.getInt(2,i)!=0 && Math.abs(pt_test.getDouble(9,i))<20 ){ //&& Math.abs(pt_test.getDouble(9,i))<20 && pt_test.getDouble(12,i)<1000
                        if(pt_test.getInt(3,i)==11) {
                          hasInstaEl_hasAIEl++;
                        }
  
                      }
                    }
  
                  }
                }
              }
            }
          }

          //has REC AI e-
          for(int j = 0; j < nrowstest; j++){
            if(pt_test.getInt(2,j)!=0 && Math.abs(pt_test.getDouble(9,j))<20 ){ //&& Math.abs(pt_test.getDouble(9,j))<20  && pt_test.getDouble(12,j)<1000
              if(pt_test.getInt(3,j)==11){
                

                double Px_o=pt_test.getDouble(4,j);
                double Py_o=pt_test.getDouble(5,j);
                double Pz_o=pt_test.getDouble(6,j);
                double P_o=Math.sqrt(Px_o*Px_o+Py_o*Py_o+Pz_o*Pz_o);
                double Theta_o = (float) Math.acos((float)Pz_o / P_o)*(180.0 / Math.PI);// Math.atan2(Math.sqrt(Px*Px+Py*Py),Pz);
                double Phi_o = (float) Math.atan2(Py_o, Px_o)*(180.0 / Math.PI);

                double var_to_check=P_o;
                if(varName=="Theta"){
                  var_to_check=Theta_o;
                } else if(varName=="Phi"){
                  var_to_check=Phi_o;
                }

                int isMatched=1;
                if(requireMatching){
                  double[] el= new double[3];
                  el[0]=Px_o;
                  el[1]=Py_o;
                  el[2]=Pz_o;
                  double[] match= new double[3];
                  isMatched=hasMatched(pt, el, match, -1);
                }

                if(isMatched!=-1){
                  if(var_to_check>=var && var_to_check<(var+step)){

                    hasAIEl++;

                    for(int i = 0; i < nrows; i++){
                      if(pt.getInt(2,i)!=0){
                        if(pt.getInt(3,i) == 11){
                          hasAIEl_hasInstaEl++;
                        }
                      }
                    }

                  }
                }

              }
            }
          }

      }

      System.out.printf("\nBinning "+varName+" [%.3f - %.3f "+varUnits+"]]\n",var,var+step);
      System.out.printf("Insta El %d & REC El %d \n", hasInstaEl,hasInstaEl_hasAIEl);
      System.out.printf("REC El %d & Insta El %d \n\n", hasAIEl,hasAIEl_hasInstaEl);

      //can get eff >1 due to more than one Insta electron per event with REC electron
      //don't really want to plot this here because we're not wanting to measure
      //increase in electron efficiency, we want to measure how many times we find
      // an electron when REC has an electron
      if(hasAIEl_hasInstaEl>hasAIEl){hasAIEl_hasInstaEl=hasAIEl;}

      StatNumber Pur = new StatNumber(hasInstaEl_hasAIEl, Math.sqrt(hasInstaEl_hasAIEl));
      StatNumber Eff = new StatNumber(hasAIEl_hasInstaEl, Math.sqrt(hasAIEl_hasInstaEl));
      StatNumber denomPur = new StatNumber(hasInstaEl, Math.sqrt(hasInstaEl));
      StatNumber denomEff = new StatNumber(hasAIEl, Math.sqrt(hasAIEl));
      Pur.divide(denomPur);
      Eff.divide(denomEff);

      gPur.addPoint(var + step / 2, Pur.number(), 0, Pur.error());
      gEff.addPoint(var + step / 2, Eff.number(), 0,Eff.error());

    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Metrics vs " + varName);
    c.draw(gEff).draw(gPur, "same");
    c.region().showLegend(0.65,0.25);
  }


  public static void main(String[] args) {

    String file="output_test_inbending_large.h5";

    plotter pter = new plotter();
    pter.plot(file,"P", "[GeV]", 1, 9.0, 1.0,true);
    //pter.plot(file, "Theta", "[Deg]", 5.0, 35, 8.,true);
    //pter.plot(file, "Phi", "[Deg]", -180, 180, 10.,true);
    //pter.plot(file,"vz", "[cm]", 1.0, 20.0, 1.0);

  }

}

