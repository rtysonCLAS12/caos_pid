import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.extratrees.networks.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;


import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.H2F;
import twig.graphics.TGCanvas;

//run with /open train_predict.java
// train_test();


void train_test(){

    int nClass=2;
    //35 when writing all HTCC mirrors out, 28 when summing the ADCs instead
    //31 with p, theta, phi
    int nVars=31; 

    //String trainDir="/Users/tyson/data_repo/trigger_data/sims/claspyth_train/for_pid/";
    String trainDir="/Users/tyson/data_repo/trigger_data/rgd/018326/for_caos_pid/";

    DataList dlt = DataList.fromCSV(trainDir+"train_fromcfpred_allNegBG2.csv",
            DataList.range(0,nVars), DataList.range(nVars,nVars+nClass));
    DataList dle = DataList.fromCSV(trainDir+"test_fromcfpred_allNegBG2.csv",
            DataList.range(0,nVars), DataList.range(nVars,nVars+1));
    DataList dlv = DataList.fromCSV(trainDir+"test_fromcfpred_allNegBG2.csv",
            DataList.range(nVars+nClass,nVars+nClass+3), DataList.range(nVars,nVars+nClass));
    dlt.shuffle();
    
    //this works nicely but unfortunately is a bit annoying to associate p to testing because of shuffling
    /*DataList dlo = DataList.fromCSV("/Users/tyson/data_repo/trigger_data/sims/claspyth_train/for_pid/train.csv",
            DataList.range(0,35), DataList.range(35,35+nClass));
    
    dlo.shuffle();
    DataList[] dls= DataList.split(dlo,0.8,0.2);
    DataList dlt=dls[0];
    DataList dle=dls[1];*/
   
    //dlt.show();
    dlt.scan();
    dlv.scan();

    ClassifierExtraTrees classifier = new ClassifierExtraTrees(nVars,nClass);

    //------ these are default values, can be commented
    classifier.setK(nVars); // set number of features in bootstrapping
    classifier.setNMin(15); // minimumu number of rows in the sample
    classifier.setNumTrees(100); // set number of trees

    classifier.train(dlt);
    classifier.export("pid_elNegBG2_fromcfpred_ert.network");

    ClassifierExtraTrees model = new ClassifierExtraTrees();
    model.load("pid_elNegBG2_fromcfpred_ert.network");


    PlotResponse(dle,model,0,"e-",nClass);
    double bestTh=findBestThreshold(dle,model,0,nClass,0.995);
    plotVarDep(dle,dlv,model,bestTh,0,2,true,0,"P","[GeV]",1,9.0,1.0);
    plotVarDep(dle,dlv,model,bestTh,0,2,true,1,"Theta","[Deg]",5.0,35.0,5.);
    plotVarDep(dle,dlv,model,bestTh,0,2,true,2,"Phi","[Deg]",-180,180,10.);

    float[] output = new float[nClass];

    /*for (int k = 0; k < 10; k++) { //dle.getList().size()
        float[] input = dle.getList().get(k).floatFirst();
        float[] desired = dle.getList().get(k).floatSecond();
        model.getOutput(input, output);
        System.out.println(Arrays.toString(input)
                + " => " + Arrays.toString(desired)
                + " => " + Arrays.toString(output));
    }*/

}

public static void PlotResponse(DataList dle,  ClassifierExtraTrees model, int elClass,String part, int nClass) {
    int NEvents = dle.getList().size();

    H1F hRespPos = new H1F(part+" in Sector", 100, 0, 1);
    hRespPos.attr().setLineColor(2);
    hRespPos.attr().setFillColor(2);
    hRespPos.attr().setLineWidth(3);
    hRespPos.attr().setTitleX("Response");
    H1F hRespNeg = new H1F("No "+part+" in Sector", 100, 0, 1);
    hRespNeg.attr().setLineColor(5);
    hRespNeg.attr().setLineWidth(3);
    hRespNeg.attr().setTitleX("Response");
    //Sort predictions into those made on the positive/or negative samples
    for(int i=0;i<NEvents;i+=1) {
        double[] input = dle.getList().get(i).getFirst();
        float[] desired = dle.getList().get(i).floatSecond();
        double output=model.evaluate(input);

        if(desired[elClass]==1) {
            hRespPos.fill(output);
        } else {
            hRespNeg.fill(output);
        }
    }

    TGCanvas c = new TGCanvas();
    
    c.setTitle("Response");
    c.draw(hRespPos).draw(hRespNeg,"same");
    c.region().showLegend(0.05, 0.95);
        
}//End of PlotResponse


//Labels col 0 is 1 if there's an e-, 0 otherwise
public static double[] getMetrics(DataList dle,  ClassifierExtraTrees model,double thresh,int elClass, int nClass){
    double[] metrics= new double[5];
    int nEvents = dle.getList().size();

    int nEls=0;
    double TP=0,FP=0,FN=0;
    for (int i = 0; i < nEvents; i++) {
        double[] input = dle.getList().get(i).getFirst();
        float[] desired = dle.getList().get(i).floatSecond();
        double output=model.evaluate(input);
        if (desired[elClass]==1) {
            nEls++;
            if (output > thresh) {
                TP++;
            } else {
                FN++;
            } 
        } else {
            if (output > thresh) {
                FP++;
            } 
        } // Check true label
    }
    double Pur=TP/(TP+FP);
    double Eff=TP/(TP+FN);
    metrics[0]=Pur;
    metrics[1]=Eff;
    metrics[2]=TP;
    metrics[3]=FP;
    metrics[4]=FN;

    /*System.out.printf("Theres %d electrons in sample\n", nEls);
    System.out.printf("L1 trigger fired %d times in sample\n", nTrig);*/
    return metrics;
}

//Labels col 0 is 1 if there's an e-, 0 otherwise
public static double[] getMetsForBin(DataList dle,DataList dlv,  ClassifierExtraTrees model,double thresh,int elClass, int nClass,int cutVar,double low,double high){
    double[] metrics = new double [2];
    int nEvents = dle.getList().size();

    double TP=0,FN=0,FP=0;
    for (int i = 0; i < nEvents; i++) {
        double[] input = dle.getList().get(i).getFirst();
        float[] vars = dlv.getList().get(i).floatFirst();
        float[] desired = dle.getList().get(i).floatSecond();
        double output=model.evaluate(input);
        if (vars[cutVar] > low && vars[cutVar]<high) {
            if (desired[elClass]==1) {
                if (output > thresh) {
                    TP++;
                } else {
                    FN++;
                } 
            } else {
                if (output > thresh) {
                    FP++;
                } 
            } // Check true label
        }
    }
    double Pur=TP/(TP+FP);
    double Eff=TP/(TP+FN);
    metrics[0]= Pur;
    metrics[1]= Eff;

    return metrics;
}

public static void plotVarDep(DataList dle,DataList dlv,  ClassifierExtraTrees model,double thresh,int elClass, int nClass,
            Boolean addPur, int cutVar, String varName, String varUnits,double low, double high,double step) {

        String yTitle="Metrics";
        if(!addPur){yTitle="Efficiency";}

        GraphErrors gEff = new GraphErrors();
        gEff.attr().setMarkerColor(2);
        gEff.attr().setMarkerSize(10);
        gEff.attr().setTitle("Level3 Efficiency");
        gEff.attr().setTitleX(varName+" "+varUnits);
        gEff.attr().setTitleY(yTitle);

        GraphErrors gPur = new GraphErrors();
        gPur.attr().setMarkerColor(5);
        gPur.attr().setMarkerSize(10);
        gPur.attr().setTitle("Level3 Purity");
        gPur.attr().setTitleX(varName+" "+varUnits);
        gPur.attr().setTitleY(yTitle);


        for (double q2=low;q2<high;q2+=step){
            double[] metrics=getMetsForBin(dle,dlv,model,thresh,elClass,nClass,cutVar,q2,q2+step);
            gPur.addPoint(q2+step/2, metrics[0], 0, 0);
            gEff.addPoint(q2+step/2, metrics[1], 0, 0);
        } // Increment threshold on response

        

        TGCanvas c = new TGCanvas();
        c.setTitle("Efficiency vs "+varName);
        c.draw(gEff);
        if(addPur){c.draw(gPur, "same");}
        c.region().axisLimitsY(gPur.getVectorY().getMin()-0.1, 1.05);
        c.region().showLegend(0.6, 0.25);
        

}

public double findBestThreshold(DataList dle,  ClassifierExtraTrees model,int elClass, int nClass,double effLow){
    
    GraphErrors gEff = new GraphErrors();
    gEff.attr().setMarkerColor(2);
    gEff.attr().setMarkerSize(10);
    gEff.attr().setTitle("Efficiency");
    gEff.attr().setTitleX("Response");
    gEff.attr().setTitleY("Metrics");
    GraphErrors gPur = new GraphErrors();
    gPur.attr().setMarkerColor(5);
    gPur.attr().setMarkerSize(10);
    gPur.attr().setTitle("Purity");
    gPur.attr().setTitleX("Response");
    gPur.attr().setTitleY("Metrics");
    double bestRespTh = 0;
    double bestPuratEffLow= 0;

    // Loop over threshold on the response
    for (double RespTh = 0.01; RespTh < 0.99; RespTh += 0.01) {
        double metrics[]=getMetrics(dle,model,RespTh,elClass,nClass);
        double Pur = metrics[0];
        double Eff = metrics[1];
        gPur.addPoint(RespTh, Pur, 0, 0);
        gEff.addPoint(RespTh, Eff, 0, 0);
        if (Eff > effLow) {
            if (Pur > bestPuratEffLow) {
                bestPuratEffLow = Pur;
                bestRespTh = RespTh;
            }
        }
    } // Increment threshold on response

    System.out.format("%n Best Purity at Efficiency above %f: %.3f at a threshold on the response of %.3f %n%n",
            effLow,bestPuratEffLow, bestRespTh);

    TGCanvas c = new TGCanvas();
    c.setTitle("Metrics vs Response");
    c.draw(gEff).draw(gPur, "same");
    c.region().showLegend(0.25, 0.25);
    //c.region().axisLimitsY(0.8, 1.01);

    return bestRespTh;
}

