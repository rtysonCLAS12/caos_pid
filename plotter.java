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

  public int[] hasMatchedAndInd(CompositeNode pt_test, double[] part, double[] matching_part, int charge){

    double Pz_o=part[2];
    int[] matched=new int[2];
    matched[0]=-1;
    matched[1]=-1;

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

  public Boolean passFid(CompositeNode cf, CompositeNode pt,int row,Boolean isREC){
    Boolean pass = true;
    for(int k = 11; k < 14; k++){
      int det =13;
      if(k<4){det=13;}
      else if (k>3 && k <7){det=14;}
      else{det=15;}
      if(isREC){det=det+1;}
      //System.out.printf("k %d & det %d\n",k,det);
      if(cf.getDouble(k,row)<9 && cf.getDouble(k,row)>0 ){ //&& pt.getDouble(det, row)>0
        pass=false;
      }
    }
    return pass;

  }


  public GraphErrors[] plot(String file,String varName, String varUnits, double low, double high, double step, Boolean requireMatching, int nEv){

    GraphErrors gEff = new GraphErrors();
    gEff.attr().setMarkerColor(2);
    gEff.attr().setMarkerSize(10);
    gEff.attr().setTitle("Efficiency");
    gEff.attr().setTitleX(varName + " " + varUnits);
    gEff.attr().setTitleY("Metrics");

    String endTitle="(Matched Tracks)";
    if(!requireMatching){endTitle="(Matched & Not Matched Tracks)";}
    

    GraphErrors gPur = new GraphErrors();
    gPur.attr().setMarkerColor(5);
    gPur.attr().setMarkerSize(10);
    gPur.attr().setTitle("Purity "+endTitle);
    gPur.attr().setTitleX(varName + " " + varUnits);
    gPur.attr().setTitleY("Metrics");


    for (double var = low; var < high; var += step) {
      //doing weird thing so i can mimick HTCC binning
      if(varName=="Theta" && var==21){step=7.0;}

      HipoReader r = new HipoReader(file);

      CompositeNode pt = new CompositeNode(32100,3,"sssifffffffsffffffffffff",1200);
      CompositeNode pt_test = new CompositeNode(32100,22,"sssifffffffsfffff",1200);
      CompositeNode cf = new CompositeNode(32100,11,"i",48);
      CompositeNode cf_test = new CompositeNode(32100,21,"i",48);
      

      Event event = new Event();
      int counter = 0;

      int hasInstaEl_hasAIEl=0, hasInstaEl=0;
      int hasAIEl=0, hasAIEl_hasInstaEl=0;


      while(r.hasNext() && counter<nEv){
          counter++;
          r.next(event);

          //System.out.println("reading");

          event.read(pt,32100,3);
          event.read(pt_test,32100,22);
          event.read(cf,32100,11);
          event.read(cf_test,32100,21);

          int nrows = pt.getRows();
          int nrowstest=pt_test.getRows();
          //int nrowstest_nAI=pt_test_nAI.getRows();

          //System.out.printf(" nrows insta %d nrows ai rec %d \n",nrows,nrowstest);

          //check if InstaRec has e-
          for(int j = 0; j < nrows; j++){

            if(pt.getInt(2,j)!=0){
              //has Insta El
              if(pt.getInt(3,j) == 11 ){//&& pt.getDouble(15,j)>0.003 && passFid(cf,j) && passFid(cf,pt,j,false)

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
                  int[] isMatchedArr=hasMatchedAndInd(pt_test, el, match, -1);
                  isMatched=isMatchedArr[0];
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
                      if(pt_test.getInt(2,i)!=0 && Math.abs(pt_test.getDouble(9,i))<12 ){ // && pt_test.getDouble(12,i)<350
                        if(pt_test.getInt(3,i)==11 ) { //&& passFid(cf_test,pt_test,i,true)
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
            if(pt_test.getInt(2,j)!=0 && Math.abs(pt_test.getDouble(9,j))<20 ){ // && pt_test.getDouble(12,j)<350
              if(pt_test.getInt(3,j)==11 ){ //&& passFid(cf_test,pt_test, j,true)
                

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
                  int[] isMatchedArr=hasMatchedAndInd(pt, el, match, -1);
                  isMatched=isMatchedArr[0];
                }

                if(isMatched!=-1){
                  if(var_to_check>=var && var_to_check<(var+step)){

                    hasAIEl++;

                    for(int i = 0; i < nrows; i++){
                      if(pt.getInt(2,i)!=0){
                        if(pt.getInt(3,i) == 11 ){ //&& passFid(cf, i) && passFid(cf,pt, i,false)
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

    GraphErrors[] gs = new GraphErrors[2];
    gs[0]=gEff;
    gs[1]=gPur;
    return gs;
  }


  public static void main(String[] args) {

    String file="output_test_inbending_large.h5"; //_newNetwork

    plotter pter = new plotter();
    GraphErrors[] gs_p_matched=pter.plot(file,"P", "[GeV]", 1, 9.0, 1.0,true,100000);//1000000
    GraphErrors[] gs_p_notmatched=pter.plot(file,"P", "[GeV]", 1, 9.0, 1.0,false,100000); //1000000
    //pter.plot(file, "Theta", "[Deg]", 5.0, 35, 8.,true);
    //pter.plot(file, "Phi", "[Deg]", -180, 180, 10.,true);
    //pter.plot(file,"vz", "[cm]", 1.0, 20.0, 1.0);

    gs_p_notmatched[1].attr().setMarkerColor(6);

    TGCanvas c = new TGCanvas();
    c.setTitle("Metrics vs P");
    c.draw(gs_p_matched[0]).draw(gs_p_matched[1], "same").draw(gs_p_notmatched[1], "same");
    c.region().showLegend(0.65,0.25);

  }

}

