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
public class physics_plotter {

  public physics_plotter() {

  }

  public double square(double a){
    return a*a;
  }

  public double getM(int pid) {
    switch (pid) {
      case 22:
        return 0;
      case 11:
        return 0.000511;
      case -11:
        return 0.000511;
      case 211:
        return 0.13957;
      case -211:
        return 0.13957;
      case 13:
        return 0.10566;
      case -13:
        return 0.10566;
      case 321:
        return 0.49368;
      case -321:
        return 0.49368;
      case 2212:
        return 0.938272;
      case 2112:
        return 0.939565;
      default:
        return -1;
    }
  }

  public String getOutString(double[] el,double[] pim,double[] pip, double[] Ms){

    double elP=Math.sqrt(square(el[0])+square(el[1])+square(el[2]));
    double elTheta = Math.acos(el[2] / elP)*180/Math.PI;// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double elPhi = Math.atan2(el[1], el[0])*180/Math.PI;

    double pimP=Math.sqrt(square(pim[0])+square(pim[1])+square(pim[2]));
    double pimTheta = Math.acos(pim[2] / pimP)*180/Math.PI;// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double pimPhi = Math.atan2(pim[1], pim[0])*180/Math.PI;

    double pipP=Math.sqrt(square(pip[0])+square(pip[1])+square(pip[2]));
    double pipTheta = Math.acos(pip[2] / pipP)*180/Math.PI;// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double pipPhi = Math.atan2(pip[1], pip[0])*180/Math.PI;

    StringBuilder csvLineBuilder = new StringBuilder();
    csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f", Ms[0],Ms[1],elP,elTheta,elPhi,pimP,pimTheta,pimPhi,pipP,pipTheta,pipPhi,Ms[2],Ms[3],Ms[4],Ms[5]));// /3
    return csvLineBuilder.toString();
  }

  public void calc_exc(double beamE,double[] el,double[] pim,double[] pip, double[] Ms){

    double elP=Math.sqrt(square(el[0])+square(el[1])+square(el[2]));
    double elE=Math.sqrt(square(elP)+square(getM(11)));

    double pimP=Math.sqrt(square(pim[0])+square(pim[1])+square(pim[2]));
    double pimE=Math.sqrt(square(pimP)+square(getM(-211)));

    double pipP=Math.sqrt(square(pip[0])+square(pip[1])+square(pip[2]));
    double pipE=Math.sqrt(square(pipP)+square(getM(211)));

    double IM = Math.sqrt(square(pimE+pipE)- ( square(pim[0]+pip[0]) + square(pim[1]+pip[1]) + square(pim[2]+pip[2]) ));

    double px_m = -1.0*(el[0]+pim[0]+pip[0]);
    double py_m = -1.0*(el[1]+pim[1]+pip[1]);
    double pz_m = beamE - (el[2]+pim[2]+pip[2]);
    double E_m = beamE + 0.938 - (elE+pipE+pimE);
    double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));
    double P_m=Math.sqrt(square(px_m)+square(py_m)+square(pz_m));
    double Theta_m = Math.acos(pz_m / P_m);// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double Phi_m = Math.atan2(py_m, px_m);
    Ms[0]=IM;
    Ms[1]=Math.sqrt(MM2);
    Ms[2]=P_m;
    Ms[3]=E_m;
    Ms[4]=Theta_m;
    Ms[5]=Phi_m;
  }

  public void calc_excK(double beamE,double[] el,double[] pim,double[] pip, double[] Ms){

    double elP=Math.sqrt(square(el[0])+square(el[1])+square(el[2]));
    double elE=Math.sqrt(square(elP)+square(getM(11)));

    double pimP=Math.sqrt(square(pim[0])+square(pim[1])+square(pim[2]));
    double pimE=Math.sqrt(square(pimP)+square(getM(-321)));

    double pipP=Math.sqrt(square(pip[0])+square(pip[1])+square(pip[2]));
    double pipE=Math.sqrt(square(pipP)+square(getM(321)));

    double IM = Math.sqrt(square(pimE+pipE)- ( square(pim[0]+pip[0]) + square(pim[1]+pip[1]) + square(pim[2]+pip[2]) ));

    double px_m = -1.0*(el[0]+pim[0]+pip[0]);
    double py_m = -1.0*(el[1]+pim[1]+pip[1]);
    double pz_m = beamE - (el[2]+pim[2]+pip[2]);
    double E_m = beamE + 0.938 - (elE+pipE+pimE);
    double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));
    double P_m=Math.sqrt(square(px_m)+square(py_m)+square(pz_m));
    double Theta_m = Math.acos(pz_m / P_m);// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double Phi_m = Math.atan2(py_m, px_m);
    Ms[0]=IM;
    Ms[1]=Math.sqrt(MM2);
    Ms[2]=P_m;
    Ms[3]=E_m;
    Ms[4]=Theta_m;
    Ms[5]=Phi_m;
  }

  public void calc_excKp(double beamE,double[] el,double[] p,double[] Km, double[] Ms){

    double elP=Math.sqrt(square(el[0])+square(el[1])+square(el[2]));
    double elE=Math.sqrt(square(elP)+square(getM(11)));

    double pP=Math.sqrt(square(p[0])+square(p[1])+square(p[2]));
    double pE=Math.sqrt(square(pP)+square(getM(2212)));

    double KmP=Math.sqrt(square(Km[0])+square(Km[1])+square(Km[2]));
    double KmE=Math.sqrt(square(KmP)+square(getM(-321)));

    double IM = Math.sqrt(square(pE+KmE)- ( square(p[0]+Km[0]) + square(p[1]+Km[1]) + square(p[2]+Km[2]) ));

    double px_m = -1.0*(el[0]+p[0]+Km[0]);
    double py_m = -1.0*(el[1]+p[1]+Km[1]);
    double pz_m = beamE - (el[2]+p[2]+Km[2]);
    double E_m = beamE + 0.938 - (elE+KmE+pE);
    double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));
    double P_m=Math.sqrt(square(px_m)+square(py_m)+square(pz_m));
    double Theta_m = Math.acos(pz_m / P_m);// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double Phi_m = Math.atan2(py_m, px_m);
    Ms[0]=IM;
    Ms[1]=Math.sqrt(MM2);
    Ms[2]=P_m;
    Ms[3]=E_m;
    Ms[4]=Theta_m;
    Ms[5]=Phi_m;
  }

  public void calc_exc_2part(double beamE,double[] el,double[] pim, double[] Ms){

    double elP=Math.sqrt(square(el[0])+square(el[1])+square(el[2]));
    double elE=Math.sqrt(square(elP)+square(getM(11)));

    double pimP=Math.sqrt(square(pim[0])+square(pim[1])+square(pim[2]));
    double pimE=Math.sqrt(square(pimP)+square(getM(-211)));

    double IM = Math.sqrt(square(pimE+elE)- ( square(pim[0]+el[0]) + square(pim[1]+el[1]) + square(pim[2]+el[2]) ));

    double px_m = -1.0*(el[0]+pim[0]);
    double py_m = -1.0*(el[1]+pim[1]);
    double pz_m = beamE - (el[2]+pim[2]);
    double E_m = beamE + 0.938 - (elE+pimE);
    double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));
    double P_m=Math.sqrt(square(px_m)+square(py_m)+square(pz_m));
    double Theta_m = Math.acos(pz_m / P_m);// Math.atan2(Math.sqrt(px*px+py*py),pz);
    double Phi_m = Math.atan2(py_m, px_m);
    Ms[0]=IM;
    Ms[1]=Math.sqrt(MM2);
    Ms[2]=P_m;
    Ms[3]=E_m;
    Ms[4]=Theta_m;
    Ms[5]=Phi_m;
  }

  public int hasMatched(CompositeNode pt_test, double[] part, double[] matching_part, int charge){

    double Pz_o=part[2];
    int matched=-1;

    double bestRes=999;
    int bestInd=-1;

    for(int i = 0; i < pt_test.getRows(); i++){
      double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
      if(pt_test.getInt(1,i)==charge && resPz<0.2){
        if(resPz<bestRes){
          bestRes=resPz;
          bestInd=i;
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

  public int[] hasMatchedAndInd(CompositeNode pt_test, double[] part, double[] matching_part, int charge){

    double Pz_o=part[2];
    int[] matched=new int[2];
    matched[0]=-1;
    matched[1]=-1;

    double bestRes=999;
    int bestInd=-1;

    for(int i = 0; i < pt_test.getRows(); i++){
      double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
      if(pt_test.getInt(1,i)==charge && resPz<0.2){
        if(resPz<bestRes){
          bestRes=resPz;
          bestInd=i;
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


  public int find_best(CompositeNode pt,int pid, double ps[],Boolean isREC){
    int nrows = pt.getRows();
    double bestP=0;
    int bestInd=-1;

    for(int i = 0; i < nrows; i++){
      int rpid=pt.getInt(3,i);

      //for instarec, any charged particle that isn't an electron is called a pion
      //do same thing for REC
      if(isREC){
        if(pt.getInt(1,i)==1){
          rpid=211;
        }
        if(pt.getInt(1,i)==-1 && rpid!=11){
          rpid=-211;
        }
      }

      if(rpid==pid){
        double Px=pt.getDouble(4,i);
        double Py=pt.getDouble(5,i);
        double Pz=pt.getDouble(6,i);
        double P=Math.sqrt(Px*Px+Py*Py+Pz*Pz);
        if(P>bestP){
          bestP=P;
          bestInd=i;
        }
      }
    }

    if(bestInd!=-1){
      ps[0]=pt.getDouble(4,bestInd);
      ps[1]=pt.getDouble(5,bestInd);
      ps[2]=pt.getDouble(6,bestInd);
    } else{
      ps[0]=0;
      ps[1]=0;
      ps[2]=0;
    }
    return bestInd;

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

  public void plot_ngamma(String file,double beamE, int nEv){

    H2F hThetaPhi = new H2F("dPhi vs dTheta",50,-1.5,1.5,50,-15,5);
    hThetaPhi.attr().setTitleX("dTheta [Degrees]");
    hThetaPhi.attr().setTitleY("dPhi [Degrees]");

    H1F hP_Insta = new H1F("Rec ID !11 & InstaRec ID 11",50,0.,10);
    hP_Insta.attr().setLineColor(2);
    hP_Insta.attr().setFillColor(2);
    hP_Insta.attr().setLineWidth(3);
    hP_Insta.attr().setTitleX("Momentum [GeV]");
    H1F hP_InstaBad = new H1F("Rec ID !11 & InstaRec ID !11",50,0.,10);
    hP_InstaBad.attr().setLineColor(3);
    hP_InstaBad.attr().setLineWidth(3);
    hP_InstaBad.attr().setTitleX("Momentum [GeV]");
    H1F hP_REC = new H1F("Rec ID !11",50,0.,10);
    hP_REC.attr().setLineColor(5);
    hP_REC.attr().setLineWidth(3);
    hP_REC.attr().setTitleX("Momentum [GeV]");

    H1F hP_SF_Insta = new H1F("Rec ID !11 & InstaRec ID 11",50,0.,10);
    hP_SF_Insta.attr().setLineColor(2);
    hP_SF_Insta.attr().setFillColor(2);
    hP_SF_Insta.attr().setLineWidth(3);
    hP_SF_Insta.attr().setTitleX("Momentum [GeV]");
    H1F hP_SF_InstaBad = new H1F("Rec ID !11 & InstaRec ID !11",50,0.,10);
    hP_SF_InstaBad.attr().setLineColor(3);
    hP_SF_InstaBad.attr().setLineWidth(3);
    hP_SF_InstaBad.attr().setTitleX("Momentum [GeV]");
    H1F hP_SF_REC = new H1F("Rec ID !11",50,0.,10);
    hP_SF_REC.attr().setLineColor(5);
    hP_SF_REC.attr().setLineWidth(3);
    hP_SF_REC.attr().setTitleX("Momentum [GeV]");

    H1F hP_HTCC_Insta = new H1F("Rec ID !11 & InstaRec ID 11",50,0.,10);
    hP_HTCC_Insta.attr().setLineColor(2);
    hP_HTCC_Insta.attr().setFillColor(2);
    hP_HTCC_Insta.attr().setLineWidth(3);
    hP_HTCC_Insta.attr().setTitleX("Momentum [GeV]");
    H1F hP_HTCC_InstaBad = new H1F("Rec ID !11 & InstaRec ID !11",50,0.,10);
    hP_HTCC_InstaBad.attr().setLineColor(3);
    hP_HTCC_InstaBad.attr().setLineWidth(3);
    hP_HTCC_InstaBad.attr().setTitleX("Momentum [GeV]");
    H1F hP_HTCC_REC = new H1F("Rec ID !11",50,0.,10);
    hP_HTCC_REC.attr().setLineColor(5);
    hP_HTCC_REC.attr().setLineWidth(3);
    hP_HTCC_REC.attr().setTitleX("Momentum [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);
    CompositeNode cf = new CompositeNode(32100,11,"i",48);
    CompositeNode cf_test = new CompositeNode(32100,21,"i",48);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext() && counter<nEv) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      event.read(cf,32100,11);
      event.read(cf_test,32100,21);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      //loop over particles, find neg particles
      //that aren't ID as electrons by REC
      for(int i = 0; i < nrows_test; i++){
        int pid_neg=pt_test.getInt(3,i);

        if(pt_test.getInt(1,i)==-1 && pid_neg!=11){
          pid_neg=-211;
        }

        if(pid_neg==-211 && passFid(cf_test, pt_test, i, true)){ //
          double Px_neg=pt_test.getDouble(4,i);
          double Py_neg=pt_test.getDouble(5,i);
          double Pz_neg=pt_test.getDouble(6,i);
          double P_neg=Math.sqrt(Px_neg*Px_neg+Py_neg*Py_neg+Pz_neg*Pz_neg);
          double Theta_neg = Math.acos(Pz_neg / P_neg)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
          double Phi_neg = Math.atan2(Py_neg, Px_neg)*(180/Math.PI);
          double[] neg = new double[3];
          neg[0]=Px_neg;
          neg[1]=Py_neg;
          neg[2]=Pz_neg;
          double nphe=pt_test.getDouble(16,i);
          double PCALE=pt_test.getDouble(13,i);

          //loop over particles, find neutral particles
          for(int j = 0; j < nrows_test; j++){
            int pid_neut=pt_test.getInt(3,j);
          
            //if(pt_test.getInt(1,j)==0){
            //  pid_neut=22;
            //}
            
            if(pid_neut==22){
              double Px_neut=pt_test.getDouble(4,j);
              double Py_neut=pt_test.getDouble(5,j);
              double Pz_neut=pt_test.getDouble(6,j);
              double P_neut=Math.sqrt(Px_neut*Px_neut+Py_neut*Py_neut+Pz_neut*Pz_neut);
              double Theta_neut = Math.acos(Pz_neut / P_neut)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
              double Phi_neut = Math.atan2(Py_neut, Px_neut)*(180/Math.PI);

              double dTheta=Theta_neut-Theta_neg;
              double dPhi=Phi_neut-Phi_neg;

              hThetaPhi.fill(dTheta, dPhi);

              double[] el=new double[3];
              int[] match=hasMatchedAndInd(pt,neg,el,-1);

              if(Math.abs(dTheta)<0.5 && match[0]!=-1 && PCALE>0.06){
                double HTCC_ADC=pt.getDouble(15,match[1]);
                if(HTCC_ADC!=0){
                  hP_REC.fill(P_neg);
                
                  if(match[0]==11){
                    hP_Insta.fill(P_neg);
                  } else{
                    hP_InstaBad.fill(P_neg);
                  }

                  if(nphe<2){

                    hP_HTCC_REC.fill(P_neg);
                
                    if(match[0]==11){
                      hP_HTCC_Insta.fill(P_neg);
                    } else{
                      hP_HTCC_InstaBad.fill(P_neg);
                    }

                  } else {

                    hP_SF_REC.fill(P_neg);
                
                    if(match[0]==11){
                      hP_SF_Insta.fill(P_neg);
                    } else{
                      hP_SF_InstaBad.fill(P_neg);
                    }

                  }

                }
              }
            }
          }
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Momentum");
    c.draw(hP_REC).draw(hP_Insta,"same").draw(hP_InstaBad, "same");
    c.region().showLegend(0.35,0.95);

    TGCanvas cHTCC = new TGCanvas();
    cHTCC.setTitle("Momentum (no HTCC)");
    cHTCC.draw(hP_HTCC_REC).draw(hP_HTCC_Insta,"same").draw(hP_HTCC_InstaBad, "same");
    cHTCC.region().showLegend(0.35,0.95);

    TGCanvas cSF = new TGCanvas();
    cSF.setTitle("Momentum (SF Problem)");
    cSF.draw(hP_SF_REC).draw(hP_SF_Insta,"same").draw(hP_SF_InstaBad, "same");
    cSF.region().showLegend(0.05,0.95);

    TGCanvas cd = new TGCanvas();
    cd.setTitle("dPhi vs dTheta");
    cd.draw(hThetaPhi);

  }


  public void plot_e2pi(String file,double beamE){

    H1F hMM_Insta = new H1F("InstaRec",50,0.5,1.5);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_REC = new H1F("REC AI",50,0.5,1.5);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("InstaRec",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("Di-pion Invariant Mass [GeV]");
    H1F hIM_REC = new H1F("REC AI",100,0,3.5);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("Di-pion Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows>2){
        double[] el=new double[3];
        double[] pim=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt,11,el,false);
        int pimInd=find_best(pt,-211,pim,false);
        int pipInd=find_best(pt,211,pip,false);
        if(elInd!=-1 && pimInd!=-1 && pipInd!=-1){
          double[] Ms= new double[6];
          calc_exc(beamE,el,pim,pip,Ms);
          hIM_Insta.fill(Ms[0]);
          if(Ms[0]<0.9 && Ms[0]>0.5){hMM_Insta.fill(Ms[1]);}
          //hMM_Insta.fill(Ms[1]);
        }
      }

      if(nrows_test>2){
        double[] el=new double[3];
        double[] pim=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt_test,11,el,true);
        int pimInd=find_best(pt_test,-211,pim,true);
        int pipInd=find_best(pt_test,211,pip,true);
        if(elInd!=-1 && pimInd!=-1 && pipInd!=-1){
          double[] Ms= new double[6];
          calc_exc(beamE,el,pim,pip,Ms);
          hIM_REC.fill(Ms[0]);
          if(Ms[0]<0.9 && Ms[0]>0.5){hMM_REC.fill(Ms[1]);}
          //hMM_REC.fill(Ms[1]);
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta).draw(hMM_REC, "same");
    c.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("Di-pion Invariant Mass");
    cIM.draw(hIM_Insta).draw(hIM_REC, "same");
    cIM.region().showLegend(0.05,0.95);
  }

  public void plot_e2pi_matchedTracks(String file,double beamE){

    H1F hMM_Insta = new H1F("All Matched e- Tracks",100,0.5,1.5);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_Insta2 = new H1F("Matched e- with different pid to REC",100,0.5,1.5);
    hMM_Insta2.attr().setLineColor(5);
    hMM_Insta2.attr().setLineWidth(3);
    hMM_Insta2.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("All Matched e- Tracks",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("Di-pion Invariant Mass [GeV]");
    H1F hIM_Insta2 = new H1F("Matched e- with different pid to REC",100,0,3.5);
    hIM_Insta2.attr().setLineColor(5);
    hIM_Insta2.attr().setLineWidth(3);
    hIM_Insta2.attr().setTitleX("Di-pion Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows>2){
        double[] el=new double[3];
        double[] pim=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt,11,el,false);
        int pimInd=find_best(pt,-211,pim,false);
        int pipInd=find_best(pt,211,pip,false);
        if(elInd!=-1 && pimInd!=-1 && pipInd!=-1){
          
          double[] matching_el=new double[3];
          double[] matching_pim=new double[3];
          double[] matching_pip=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);
          int pipMatch=hasMatched(pt_test,pip,matching_pip,1);
          int pimMatch=hasMatched(pt_test,pim,matching_pim,-1);
          if(elMatch!=-1 && pimMatch!=-1 && pipMatch!=-1){ //&& pimMatch!=-1 && pipMatch!=-1
            double[] Ms= new double[6];
            calc_exc(beamE,matching_el,matching_pim,matching_pip,Ms); //matching_
            hIM_Insta.fill(Ms[0]);
            if(Ms[2]<2){hMM_Insta.fill(Ms[1]);}
            
            if(elMatch!=11){
              hIM_Insta2.fill(Ms[0]);
              if(Ms[2]<2){hMM_Insta2.fill(Ms[1]);}
            }

          }
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta).draw(hMM_Insta2,"same");
    c.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("Di-pion Invariant Mass");
    cIM.draw(hIM_Insta).draw(hIM_Insta2,"same");
    cIM.region().showLegend(0.45,0.95);
  }

  public void plot_epip_matchedTracks(String file,double beamE){

    H1F hMM_Insta = new H1F("All Matched e- Tracks",50,0.5,1.5);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_Insta2 = new H1F("Matched e- with different pid to REC",50,0.5,1.5);
    hMM_Insta2.attr().setLineColor(5);
    hMM_Insta2.attr().setLineWidth(3);
    hMM_Insta2.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("All Matched e- Tracks",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("e' pi+ Invariant Mass [GeV]");
    H1F hIM_Insta2 = new H1F("Matched e- with different pid to REC",100,0,3.5);
    hIM_Insta2.attr().setLineColor(5);
    hIM_Insta2.attr().setLineWidth(3);
    hIM_Insta2.attr().setTitleX("e' pi+ Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows>1){
        double[] el=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt,11,el,false);
        int pipInd=find_best(pt,211,pip,false);
        if(elInd!=-1 && pipInd!=-1){
          
          double[] matching_el=new double[3];
          double[] matching_pip=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);
          int pipMatch=hasMatched(pt_test,pip,matching_pip,1);
          if(elMatch!=-1 && pipMatch!=-1){ //&& pipMatch!=-1
            double[] Ms= new double[6];
            calc_exc_2part(beamE,matching_el,matching_pip,Ms);//matching_
            hIM_Insta.fill(Ms[0]);
            if(Ms[2]<1.){hMM_Insta.fill(Ms[1]);}
            if(elMatch!=11){
              hIM_Insta2.fill(Ms[0]);
              if(Ms[2]<1.){hMM_Insta2.fill(Ms[1]);}
            }

          }
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta).draw(hMM_Insta2,"same");
    c.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("e' pi+ Invariant Mass");
    cIM.draw(hIM_Insta).draw(hIM_Insta2,"same");
    cIM.region().showLegend(0.45,0.95);
  }

  public void plot_epip(String file,double beamE){

    H1F hMM_Insta = new H1F("InstaRec",50,0.5,1.5);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_REC = new H1F("REC AI",50,0.5,1.5);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("InstaRec",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("e' pi+ Invariant Mass [GeV]");
    H1F hIM_REC = new H1F("REC AI",100,0,3.5);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("e' pi+ Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows>1){
        double[] el=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt,11,el,false);
        int pipInd=find_best(pt,211,pip,false);
        if(elInd!=-1 && pipInd!=-1){
          double[] Ms= new double[6];
          calc_exc_2part(beamE,el,pip,Ms);
          hMM_Insta.fill(Ms[1]);
          hIM_Insta.fill(Ms[0]);
        }
      }

      
      if(nrows_test>1){
        double[] el=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt_test,11,el,true);
        int pipInd=find_best(pt_test,211,pip,true);
        if(elInd!=-1 && pipInd!=-1){
          double[] Ms= new double[6];
          calc_exc_2part(beamE,el,pip,Ms);
          hMM_REC.fill(Ms[1]);
          hIM_REC.fill(Ms[0]);
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta).draw(hMM_REC, "same");
    c.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("e' pi+ Invariant Mass");
    cIM.draw(hIM_Insta).draw(hIM_REC, "same");
    cIM.region().showLegend(0.05,0.95);

  }

  public void plot_epipREC(String file,double beamE){

    H1F hMM_REC = new H1F("REC AI",50,0.5,1.5);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_REC = new H1F("REC AI",100,0,3.5);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("e' pi+ Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();
      
      if(nrows_test>1){
        double[] el=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt_test,11,el,false);
        int pipInd=find_best(pt_test,211,pip,false);
        if(elInd!=-1 && pipInd!=-1){
          double[] Ms= new double[6];
          calc_exc_2part(beamE,el,pip,Ms);
          //if(Ms[2]<1.0){hMM_REC.fill(Ms[1]);}
          hMM_REC.fill(Ms[1]);
          hIM_REC.fill(Ms[0]);
          
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_REC);
    c.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("e' pi+ Invariant Mass");
    cIM.draw(hIM_REC);
    cIM.region().showLegend(0.05,0.95);

  }

  //best to get estimate of purity of pid only
  public void plot_e2piREC(String file,double beamE){

    H1F hMM_Insta = new H1F("REC AI e' pi- pi+ InstaRec e' e- pi+",100,0.5,3.5);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hnoMM_Insta = new H1F("REC AI & InstaRec e' pi- pi+",100,0.5,3.5);
    hnoMM_Insta.attr().setLineColor(4);
    hnoMM_Insta.attr().setLineWidth(3);
    hnoMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_REC = new H1F("REC AI e' pi- pi+",100,0.5,3.5);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("REC AI e' pi- pi+ InstaRec e' e- pi+",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("Di-pion Invariant Mass [GeV]");
    H1F hnoIM_Insta = new H1F("REC AI & InstaRec e' pi- pi+",100,0,3.5);
    hnoIM_Insta.attr().setLineColor(4);
    hnoIM_Insta.attr().setLineWidth(3);
    hnoIM_Insta.attr().setTitleX("Di-pion Invariant Mass [GeV]");
    H1F hIM_REC = new H1F("REC AI e' pi- pi+",100,0,3.5);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("Di-pion Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows_test>2){
        double[] el=new double[3];
        double[] pim=new double[3];
        double[] pip=new double[3];
        int elInd=find_best(pt_test,11,el,false);
        int pimInd=find_best(pt_test,-211,pim,false);
        int pipInd=find_best(pt_test,211,pip,false);
        if(elInd!=-1 && pimInd!=-1 && pipInd!=-1){
          double[] matching_pim=new double[3];
          int pimMatch=hasMatched(pt,pim,matching_pim,-1);
          double[] Ms= new double[6];
          calc_exc(beamE,el,pim,pip,Ms);
          hIM_REC.fill(Ms[0]);
          //if(Ms[0]<0.9 && Ms[0]>0.5){hMM_REC.fill(Ms[1]);}
          hMM_REC.fill(Ms[1]);
          if(pimMatch==11){
            hIM_Insta.fill(Ms[0]);
            //if(Ms[0]<0.9 && Ms[0]>0.5){hMM_Insta.fill(Ms[1]);}
            hMM_Insta.fill(Ms[1]);

          } else{
            hnoIM_Insta.fill(Ms[0]);
            //if(Ms[0]<0.9 && Ms[0]>0.5){hnoMM_Insta.fill(Ms[1]);}
            hnoMM_Insta.fill(Ms[1]);
          }
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta);
    c.region().showLegend(0.05,0.95);

    TGCanvas c2 = new TGCanvas();
    c2.setTitle("Missing Mass");
    c2.draw(hMM_REC).draw(hnoMM_Insta,"same");
    c2.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("Di-pion Invariant Mass");
    cIM.draw(hIM_Insta);
    cIM.region().showLegend(0.05,0.95);

    TGCanvas cIM2 = new TGCanvas();
    cIM2.setTitle("Di-pion Invariant Mass");
    cIM2.draw(hIM_REC).draw(hnoIM_Insta,"same");
    cIM2.region().showLegend(0.05,0.95);
  }

  //best to get estimate of purity of pid and tracking
  public void plot_2piREC(String file,double beamE){

    H1F hMM_Insta = new H1F("REC AI pi+ pi- InstaRec e' ",50,0.5,2.0);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_Insta_wCut = new H1F("REC AI pi+ pi- InstaRec e' (Cut on IM)",50,0.5,2.0);
    hMM_Insta_wCut.attr().setLineColor(6);
    hMM_Insta_wCut.attr().setLineWidth(3);
    hMM_Insta_wCut.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("REC AI pi+ pi- InstaRec e'",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("Di-pion Invariant Mass [GeV]");
    H1F hIM_REC = new H1F("REC AI pi+ pi-",100,0,3.5);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("Di-pion Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    Path filePath = Paths.get("forplot.csv");

    // Delete the file if it already exists
    try {
      Files.deleteIfExists(filePath);
    } catch (IOException e) {
      e.printStackTrace();
    }

    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath))) {

      while (r.hasNext()) {
        counter++;
        r.next(event);

        // System.out.println("reading");

        event.read(pt, 32100, 3);
        event.read(pt_test, 32100, 22);
        int nrows = pt.getRows();
        int nrows_test = pt_test.getRows();

        if(nrows>2 && nrows_test>1){
          
          double[] pim=new double[3];
          double[] pip=new double[3];
          double[] el=new double[3];
          int pimInd=find_best(pt,-211,pim,false);
          int pipInd=find_best(pt,211,pip,false);
          int elInd=find_best(pt,11,el,false);
          double[] REC_pim=new double[3];
          double[] REC_pip=new double[3];
          double[] REC_el=new double[3];
          int REC_pimInd=find_best(pt_test,-211,REC_pim,false);
          int REC_pipInd=find_best(pt_test,211,REC_pip,false);
          int REC_elInd=find_best(pt_test,11,REC_el,false);
          if( elInd!=-1 && REC_pimInd!=-1 && REC_pipInd!=-1 && REC_elInd==-1){ //pimInd!=-1 && pipInd!=-1 &&
          
            double[] matching_el=new double[3];
            int elMatch=hasMatched(pt_test,el,matching_el,-1);

            double[] Ms= new double[6];
            calc_exc(beamE,el,REC_pim,REC_pip,Ms);
            hIM_REC.fill(Ms[0]);
            if(elMatch==-1){
              hIM_Insta.fill(Ms[0]);
              hMM_Insta.fill(Ms[1]);
              if(Ms[0]>0.5 && Ms[0]<1.1){hMM_Insta_wCut.fill(Ms[1]);}

              writer.println(getOutString(el,pim,pip,Ms));
              
            }
          }
        }
      
      }
    }catch (IOException e) {
      e.printStackTrace();
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta).draw(hMM_Insta_wCut,"same");
    c.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("Di-pion Invariant Mass");
    cIM.draw(hIM_Insta).draw(hIM_REC,"same");
    cIM.region().showLegend(0.05,0.95);

    TGCanvas cIM2 = new TGCanvas();
    cIM2.setTitle("Di-pion Invariant Mass");
    cIM2.draw(hIM_Insta);
    cIM2.region().showLegend(0.05,0.95);
  }

  //best to get estimate of purity of pid and tracking
  public void plot_2piREC_noWrite(String file,double beamE){

    H1F hMM_Insta = new H1F("REC AI pi+ pi- InstaRec unmatched e' ",50,0.5,2.0);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_Insta_wCut = new H1F("REC AI pi+ pi- InstaRec  unmatched e' (Cut on IM)",50,0.5,2.0);
    hMM_Insta_wCut.attr().setLineColor(6);
    hMM_Insta_wCut.attr().setLineWidth(3);
    hMM_Insta_wCut.attr().setTitleX("Missing Mass [GeV]");

    H1F hMM_REC = new H1F("REC AI pi+ pi- No InstaRec unmatched e' ",50,0.5,2.0);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");
    H1F hMM_REC_wCut = new H1F("REC AI pi+ pi- No InstaRec unmatched e' (Cut on IM)",50,0.5,2.0);
    hMM_REC_wCut.attr().setLineColor(3);
    hMM_REC_wCut.attr().setLineWidth(3);
    hMM_REC_wCut.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("REC AI pi+ pi- InstaRec unmatched e'",100,0,3.5);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("Di-pion Invariant Mass [GeV]");
    H1F hIM_REC = new H1F("REC AI pi+ pi- No InstaRec unmatched e'",100,0,3.5);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("Di-pion Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      /*if(nrows_test>1){
        
        double[] pim=new double[3];
        double[] pip=new double[3];
        int pimInd=find_best(pt_test,-211,pim,false);
        int pipInd=find_best(pt_test,211,pip,false);
        if(pimInd!=-1 && pipInd!=-1){
          double[] Insta_el=new double[3];
          int elInstaInd=find_best(pt,11,Insta_el,false);
          double[] matching_Insta_el=new double[3];
          int elInstaMatch=hasMatched(pt_test,Insta_el,matching_Insta_el,-1);

          double[] Ms= new double[6];
          calc_exc(beamE,Insta_el,pim,pip,Ms);

          hIM_REC.fill(Ms[0]);
          if(elInstaInd!=-1 && elInstaMatch==-1){
            hIM_Insta.fill(Ms[0]);
            hMM_Insta.fill(Ms[1]);
          }
        }
      }*/

     /*if(nrows>2 && nrows_test>1){
        
        double[] pim=new double[3];
        double[] pip=new double[3];
        double[] el=new double[3];
        int pimInd=find_best(pt,-211,pim,false);
        int pipInd=find_best(pt,211,pip,false);
        int elInd=find_best(pt,11,el,false);
        if(pimInd!=-1 && pipInd!=-1 && elInd!=-1){
         
          double[] matching_el=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);
          double[] matching_pim=new double[3];
          int pimMatch=hasMatched(pt_test,pim,matching_pim,-1);
          double[] matching_pip=new double[3];
          int pipMatch=hasMatched(pt_test,pip,matching_pip,1);

          if(pimMatch!=-1 && pipMatch!=-1){
            double[] Ms= new double[6];
            calc_exc(beamE,el,matching_pim,matching_pip,Ms);

            hIM_REC.fill(Ms[0]);
            if(elMatch==-1){
              hIM_Insta.fill(Ms[0]);
              hMM_Insta.fill(Ms[1]);
              if(Ms[0]>0.5 && Ms[0]<1.1){hMM_Insta_wCut.fill(Ms[1]);}
              
            }

          }
          
        }
      }*/

      if(nrows>0 && nrows_test>1){
        
        double[] pim=new double[3];
        double[] pip=new double[3];
        double[] el=new double[3];
        int pimInd=find_best(pt,-211,pim,false);
        int pipInd=find_best(pt,211,pip,false);
        int elInd=find_best(pt,11,el,false);
        double[] REC_pim=new double[3];
        double[] REC_pip=new double[3];
        double[] REC_el=new double[3];
        int REC_pimInd=find_best(pt_test,-211,REC_pim,false);
        int REC_pipInd=find_best(pt_test,211,REC_pip,false);
        int REC_elInd=find_best(pt_test,11,REC_el,false);
        if( REC_pimInd!=-1 && REC_pipInd!=-1 ){ //pimInd!=-1 && pipInd!=-1 && REC_elInd==-1
         
          double[] matching_el=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);

          double[] Ms= new double[6];
          calc_exc(beamE,el,REC_pim,REC_pip,Ms);
          if(elInd!=-1 && elMatch==-1){
            hIM_Insta.fill(Ms[0]);
            hMM_Insta.fill(Ms[1]);
            if(Ms[0]>0.6 && Ms[0]<0.9){hMM_Insta_wCut.fill(Ms[1]);}
          } else{
            hIM_REC.fill(Ms[0]);
            hMM_REC.fill(Ms[1]);
            if(Ms[0]>0.6 && Ms[0]<0.9){hMM_REC_wCut.fill(Ms[1]);}
          }

          
          
        }
      }

      
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta).draw(hMM_Insta_wCut,"same");
    c.region().showLegend(0.05,0.95);

    TGCanvas c2 = new TGCanvas();
    c2.setTitle("Missing Mass");
    c2.draw(hMM_REC).draw(hMM_REC_wCut,"same");
    c2.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("Di-pion Invariant Mass");
    cIM.draw(hIM_Insta);
    cIM.region().showLegend(0.05,0.95);

    TGCanvas cIM2 = new TGCanvas();
    cIM2.setTitle("Di-pion Invariant Mass");
    cIM2.draw(hIM_REC);
    cIM2.region().showLegend(0.05,0.95);
  }

   //best to get estimate of purity of pid and tracking
   public void plot_piREC(String file,double beamE){

    H1F hMM_Insta = new H1F("REC AI pi+ InstaRec e' pi+",50,0.5,2.0);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");


    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows>1 && nrows_test>0){
        
        double[] pip=new double[3];
        double[] el=new double[3];
        int pipInd=find_best(pt,211,pip,false);
        int elInd=find_best(pt,11,el,false);
        double[] REC_pip=new double[3];
        double[] REC_el=new double[3];
        int REC_pipInd=find_best(pt_test,211,REC_pip,false);
        int REC_elInd=find_best(pt_test,11,REC_el,false);
        if(pipInd!=-1 && elInd!=-1 && REC_pipInd!=-1 && REC_elInd==-1){
         
          double[] matching_el=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);

          double[] Ms= new double[6];
          calc_exc_2part(beamE, el, REC_pip, Ms);
          if(elMatch==-1){
            hMM_Insta.fill(Ms[1]);
          }
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta);
    c.region().showLegend(0.05,0.95);
  }

  //best to get estimate of purity of pid and tracking
  public void plot_2KREC(String file,double beamE){

    H1F hMM_Insta = new H1F("REC AI K+ K- InstaRec unmatched e'",100,0.5,3.5);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");

    H1F hMM_REC = new H1F("REC AI K+ K- No InstaRec unmatched e'",100,0.5,3.5);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("REC AI K+ K- InstaRec unmatched e'",30,0.9,1.2);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("K+ K- Invariant Mass [GeV]");

    H1F hIM_REC = new H1F("REC AI K+ K- No InstaRec unmatched e'",30,0.9,1.2);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("K+ K- Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      /*if(nrows_test>1){
        
        double[] pim=new double[3];
        double[] pip=new double[3];
        int pimInd=find_best(pt_test,-321,pim,false);
        int pipInd=find_best(pt_test,321,pip,false);
        if(pimInd!=-1 && pipInd!=-1){
          double[] Insta_el=new double[3];
          int elInstaInd=find_best(pt,11,Insta_el,false);
          double[] matching_Insta_el=new double[3];
          int elInstaMatch=hasMatched(pt_test,Insta_el,matching_Insta_el,-1);

          double[] Ms= new double[6];
          calc_excK(beamE,Insta_el,pim,pip,Ms);

          hIM_REC.fill(Ms[0]);
          if(elInstaInd!=-1 && elInstaMatch==-1){
            hIM_Insta.fill(Ms[0]);
            hMM_Insta.fill(Ms[1]);
          }
        }
      }*/

      if(nrows>0 && nrows_test>1){
        
        double[] pim=new double[3];
        double[] pip=new double[3];
        double[] el=new double[3];
        int pimInd=find_best(pt,-321,pim,false);
        int pipInd=find_best(pt,321,pip,false);
        int elInd=find_best(pt,11,el,false);
        double[] REC_pim=new double[3];
        double[] REC_pip=new double[3];
        double[] REC_el=new double[3];
        int REC_pimInd=find_best(pt_test,-321,REC_pim,false);
        int REC_pipInd=find_best(pt_test,321,REC_pip,false);
        int REC_elInd=find_best(pt_test,11,REC_el,false);
        if(REC_pimInd!=-1 && REC_pipInd!=-1){ //&& REC_elInd==-1
         
          double[] matching_el=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);

          double[] Ms= new double[6];
          calc_excK(beamE,el,REC_pim,REC_pip,Ms);
          
          if(elInd!=-1 && elMatch==-1){
            hIM_Insta.fill(Ms[0]);
            hMM_Insta.fill(Ms[1]);
      
          } else{
            hIM_REC.fill(Ms[0]);
            hMM_REC.fill(Ms[1]);
          }

          
          
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta);
    c.region().showLegend(0.05,0.95);

    TGCanvas c2 = new TGCanvas();
    c2.setTitle("Missing Mass");
    c2.draw(hMM_REC);
    c2.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("K+ K- Invariant Mass");
    cIM.draw(hIM_Insta);
    cIM.region().showLegend(0.05,0.95);

    TGCanvas cIM2 = new TGCanvas();
    cIM2.setTitle("K+ K- Invariant Mass");
    cIM2.draw(hIM_REC);
    cIM2.region().showLegend(0.05,0.95);
  }

  //best to get estimate of purity of pid and tracking
  public void plot_pKREC(String file,double beamE){

    H1F hMM_Insta = new H1F("REC AI p K- InstaRec unmatched e'",100,0,1);
    hMM_Insta.attr().setLineColor(2);
    hMM_Insta.attr().setFillColor(2);
    hMM_Insta.attr().setLineWidth(3);
    hMM_Insta.attr().setTitleX("Missing Mass [GeV]");

    H1F hMM_REC = new H1F("REC AI p K- No InstaRec unmatched e'",100,0,1);
    hMM_REC.attr().setLineColor(5);
    hMM_REC.attr().setLineWidth(3);
    hMM_REC.attr().setTitleX("Missing Mass [GeV]");

    H1F hIM_Insta = new H1F("REC AI p K- InstaRec unmatched e'",30,1.4,1.7);
    hIM_Insta.attr().setLineColor(2);
    hIM_Insta.attr().setFillColor(2);
    hIM_Insta.attr().setLineWidth(3);
    hIM_Insta.attr().setTitleX("p K- Invariant Mass [GeV]");

    H1F hIM_REC = new H1F("REC AI p K- No InstaRec unmatched e'",30,1.4,1.7);
    hIM_REC.attr().setLineColor(5);
    hIM_REC.attr().setLineWidth(3);
    hIM_REC.attr().setTitleX("p K- Invariant Mass [GeV]");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100, 3, "sssifffffffsffffffffffff", 1200);
    CompositeNode pt_test = new CompositeNode(32100, 22, "sssifffffffsfffff", 1200);

    Event event = new Event();
    int counter = 0;

    int hasInstaEl_hasAIEl = 0, hasInstaEl = 0;
    int hasAIEl = 0, hasAIEl_hasInstaEl = 0;

    while (r.hasNext()) {
      counter++;
      r.next(event);

      // System.out.println("reading");

      event.read(pt, 32100, 3);
      event.read(pt_test, 32100, 22);
      int nrows = pt.getRows();
      int nrows_test = pt_test.getRows();

      if(nrows>0 && nrows_test>1){
        
        double[] el=new double[3];
        int elInd=find_best(pt,11,el,false);
        double[] REC_p=new double[3];
        double[] REC_Km=new double[3];
        int REC_pInd=find_best(pt_test,2212,REC_p,false);
        int REC_KmInd=find_best(pt_test,-321,REC_Km,false);
        if(REC_KmInd!=-1 && REC_pInd!=-1 ){ //&& REC_elInd==-1
         
          double[] matching_el=new double[3];
          int elMatch=hasMatched(pt_test,el,matching_el,-1);

          double[] Ms= new double[6];
          calc_excKp(beamE,el,REC_p,REC_Km,Ms);
          
          if(elInd!=-1 && elMatch==-1){
            hIM_Insta.fill(Ms[0]);
            hMM_Insta.fill(Ms[1]);
      
          } else{
            hIM_REC.fill(Ms[0]);
            hMM_REC.fill(Ms[1]);
          }
          
        }
      }
     
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Missing Mass");
    c.draw(hMM_Insta);
    c.region().showLegend(0.05,0.95);

    TGCanvas c2 = new TGCanvas();
    c2.setTitle("Missing Mass");
    c2.draw(hMM_REC);
    c2.region().showLegend(0.05,0.95);

    TGCanvas cIM = new TGCanvas();
    cIM.setTitle("p K- Invariant Mass");
    cIM.draw(hIM_Insta);
    cIM.region().showLegend(0.05,0.95);

    TGCanvas cIM2 = new TGCanvas();
    cIM2.setTitle("p K- Invariant Mass");
    cIM2.draw(hIM_REC);
    cIM2.region().showLegend(0.05,0.95);
  }


  public static void main(String[] args) {

    String file="output_test_inbending_large.h5";

    physics_plotter pter = new physics_plotter();
    //pter.plot_e2pi(file,10.547);
    //pter.plot_epip(file,10.547);

    
    //pter.plot_e2pi_matchedTracks(file,10.547);
    //pter.plot_epip_matchedTracks(file,10.547);

    //pter.plot_epipREC(file,10.547);

    //best to get estimate of purity of pid only
    //pter.plot_e2piREC(file,10.547);
    pter.plot_ngamma(file, 10.547,1000000);

    
    //pter.plot_2piREC_noWrite(file,10.547);
    //ter.plot_2KREC(file,10.547);
   //pter.plot_piREC(file,10.547);
   //pter.plot_pKREC(file,10.547);

  }

}

