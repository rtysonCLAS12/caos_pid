import j4ml.data.*;
import j4ml.deepnetts.*;
import deepnetts.net.layers.activation.*;

import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;


import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.H2F;
import twig.graphics.TGCanvas;

//run with /open CF_train_predict.java
// train_test();


void train_test_cf(){

    int nClass=9;

    //String trainDir="/Users/tyson/data_repo/trigger_data/sims/claspyth_train/for_pid/";
    String trainDir="/Users/tyson/data_repo/trigger_data/rgd/018326/for_caos_pid/";

    DataList dlt = DataList.fromCSV(trainDir+"train_cf.csv",
            DataList.range(0,6), DataList.range(6,6+nClass));
    DataList dle = DataList.fromCSV(trainDir+"test_cf.csv",
            DataList.range(0,6), DataList.range(6,6+nClass));
    DataList dlv = DataList.fromCSV(trainDir+"test_cf.csv",
    DataList.range(6+nClass,6+nClass+3), DataList.range(6,6+nClass));
    dlt.shuffle();
    
    dlt.scan();

    DeepNettsNetwork regression = new DeepNettsNetwork();
    regression.activation(ActivationType.TANH); // or ActivationType.TANH
    regression.outputActivation(ActivationType.LINEAR);
    regression.init(new int[]{6,12,36,27,18,9});
    regression.train(dlt,1000);

    regression.save("cf_el.network");

    EJMLModel model = new EJMLModel("cf_el.network", ModelType.TANH_LINEAR);

    System.out.println("network structure: " + model.summary());
    System.out.println("\n\nRunning Inference:\n------------");

    String[] dets=new String[]{"PCAL","ECIN","ECOUT"};

    for(String det:dets){
        
        calcMSE(dle, model,det);
        //plotCFExamples(dle,model,det,1);
        plotConvertedLs(dle,det);

        /*plotCFVarDep(dle,dlv,model,det,0,"P","[GeV]",1,9.0,1.0);
        plotCFVarDep(dle,dlv,model,det,1,"Theta","[Deg]",5.0,35.0,5.);
        plotCFVarDep(dle,dlv,model,det,2,"Phi","[Deg]",-180,180,10.);*/
    }

}

public void plotConvertedLs(DataList dle, String dets){


    int nEvents = dle.getList().size();

    int out_start=0;
    double end=68.5;
    int n=72;
    if(dets=="ECIN"){
        out_start=3;
        end=36.5;
        n=40;
    }
    if(dets=="ECOUT"){
        out_start=6;
        end=36.5;
        n=40;
    }

    H1F hStripsU = new H1F("U Views", n,-3.5,end);
    hStripsU.attr().setLineColor(2);
    hStripsU.attr().setLineWidth(3);
    hStripsU.attr().setTitleX(dets+" Strips");

    H1F hStripsV = new H1F("V Views", n,-3.5,end);
    hStripsV.attr().setLineColor(5);
    hStripsV.attr().setLineWidth(3);
    hStripsV.attr().setTitleX(dets+" Strips");

    H1F hStripsW = new H1F("W Views", n,-3.5,end);
    hStripsW.attr().setLineColor(3);
    hStripsW.attr().setLineWidth(3);
    hStripsW.attr().setTitleX(dets+" Strips");

    H1F hLsU = new H1F("U Views",110,-10,500);
    hLsU.attr().setLineColor(2);
    hLsU.attr().setLineWidth(3);
    hLsU.attr().setTitleX(dets+" Ls");

    H1F hLsV = new H1F("V Views",110,-10,500);
    hLsV.attr().setLineColor(5);
    hLsV.attr().setLineWidth(3);
    hLsV.attr().setTitleX(dets+" Ls");

    H1F hLsW = new H1F("W Views",110,-10,500);
    hLsW.attr().setLineColor(3);
    hLsW.attr().setLineWidth(3);
    hLsW.attr().setTitleX(dets+" Ls");

    for(int i=0;i<nEvents;i+=1) {
        float[] desired = dle.getList().get(i).floatSecond();
        for(int j=0;j<9;j++){
            desired[j]=desired[j]*500;
        }
        hLsU.fill(desired[out_start+0]);
        hLsV.fill(desired[out_start+1]);
        hLsW.fill(desired[out_start+2]);
        int[] desired_strips=convLtoStrip(desired);
        hStripsU.fill(desired_strips[out_start+0]);
        hStripsV.fill(desired_strips[out_start+1]);
        hStripsW.fill(desired_strips[out_start+2]);
    }

    TGCanvas c = new TGCanvas();
    c.setTitle(dets+" Strips");
    c.draw(hStripsU).draw(hStripsV,"same").draw(hStripsW,"same");
    c.region().showLegend(0.05, 0.95);

    TGCanvas cl = new TGCanvas();
    cl.setTitle(dets+" Ls");
    cl.draw(hLsU).draw(hLsV,"same").draw(hLsW,"same");
    cl.region().showLegend(0.05, 0.95);

}

public void plotCFExamples(DataList dle,  EJMLModel model, String dets, int nExamples){

    int out_start=0;
    int n=68;
    double end=n*3+0.5;
    if(dets=="ECIN"){
        out_start=3;
        end=n*3+0.5;
        n=36;
    }
    if(dets=="ECOUT"){
        out_start=6;
        end=n*3+0.5;
        n=36;
    }

    float[] output_l = new float[9];

    for(int i=0;i<nExamples;i+=1) {

        H1F hTrue = new H1F("True", n*3+1, -0.5, end);
        hTrue.attr().setLineColor(2);
        hTrue.attr().setFillColor(2);
        hTrue.attr().setLineWidth(3);
        hTrue.attr().setTitleX(dets+" Cluster Local Position");

        H1F hPred = new H1F("Pred", n*3+1,-0.5,end);
        hPred.attr().setLineColor(5);
        hPred.attr().setLineWidth(3);
        hPred.attr().setTitleX(dets+" Cluster Local Position");

        float[] input = dle.getList().get(i).floatFirst();
        float[] desired_l = dle.getList().get(i).floatSecond();
        model.getOutput(input, output_l);
        for(int j=0;j<9;j++){
            desired_l[j]=desired_l[j]*500;
            output_l[j]=output_l[j]*500;
        }
        int[] desired=convLtoStrip(desired_l);
        int[] output=convLtoStrip(output_l);

        hTrue.fill(desired[out_start+0]);
        hTrue.fill(desired[out_start+1]+n);
        hTrue.fill(desired[out_start+2]+n*2);

        hPred.fill(output[out_start+0]);
        hPred.fill(output[out_start+1]+n);
        hPred.fill(output[out_start+2]+n*2);
        TGCanvas c = new TGCanvas();
        c.setTitle(dets+" Cluster Position");
        c.draw(hTrue).draw(hPred,"same");
        c.region().showLegend(0.05, 0.95);
        

    }

}

public void calcMSE(DataList dle,  EJMLModel model, String dets){
    int nEvents = dle.getList().size();
    float[] output_l = new float[9];

    int out_start=0;
    double end=68.5;
    int n=69;
    if(dets=="ECIN"){
        out_start=3;
        end=36.5;
        n=37;
    }
    if(dets=="ECOUT"){
        out_start=6;
        end=36.5;
        n=37;
    }

    float adifU=0,adifV=0,adifW=0;

    H2F hDifTrueU = new H2F("U View", n,-0.5,end,11,-5.5,5.5);
    hDifTrueU.attr().setTitleY("U View Position Difference [strips]");
    hDifTrueU.attr().setTitleX("U View  True Position");

    H2F hDifTrueV = new H2F("V View", n,-0.5,end,11,-5.5,5.5);
    hDifTrueV.attr().setTitleY("V View Position Difference [strips]");
    hDifTrueV.attr().setTitleX("V View  True Position");

    H2F hDifTrueW = new H2F("W View", n,-0.5,end,11,-5.5,5.5);
    hDifTrueW.attr().setTitleY("W View Position Difference [strips]");
    hDifTrueW.attr().setTitleX("W View  True Position");

    H1F hDifU = new H1F("U View", 11,-5.5,5.5);
    hDifU.attr().setLineColor(2);
    hDifU.attr().setLineWidth(3);
    hDifU.attr().setTitleX("Position Difference [strips]");

    H1F hDifV = new H1F("V View", 11,-5.5,5.5);
    hDifV.attr().setLineColor(5);
    hDifV.attr().setLineWidth(3);
    hDifV.attr().setTitleX("Position Difference [strips]");

    H1F hDifW = new H1F("W View", 11,-5.5,5.5);
    hDifW.attr().setLineColor(3);
    hDifW.attr().setLineWidth(3);
    hDifW.attr().setTitleX("Position Difference [strips]");

    for(int i=0;i<nEvents;i+=1) {
        float[] input = dle.getList().get(i).floatFirst();
        float[] desired_l = dle.getList().get(i).floatSecond();
        model.getOutput(input, output_l);
        for(int j=0;j<9;j++){
            desired_l[j]=desired_l[j]*500;
            output_l[j]=output_l[j]*500;
        }
        int[] desired=convLtoStrip(desired_l);
        int[] output=convLtoStrip(output_l);

        adifU+=(desired[out_start+0]-output[out_start+0]);
        adifV+=(desired[out_start+1]-output[out_start+1]);
        adifW+=(desired[out_start+2]-output[out_start+2]);
        hDifU.fill((desired[out_start+0]-output[out_start+0]));
        hDifV.fill((desired[out_start+1]-output[out_start+1]));
        hDifW.fill((desired[out_start+2]-output[out_start+2]));
        hDifTrueU.fill(desired[out_start+0],(desired[out_start+0]-output[out_start+0]));
        hDifTrueV.fill(desired[out_start+1],(desired[out_start+1]-output[out_start+1]));
        hDifTrueW.fill(desired[out_start+2],(desired[out_start+2]-output[out_start+2]));

        /*if((desired[0]-output[0])>0.2 && i<100){
            System.out.printf("el nb %d U %f V %f W %f \n",i,desired[0],desired[1],desired[2]);
        }*/
    }
    adifU=adifU/nEvents;
    adifV=adifV/nEvents;
    adifW=adifW/nEvents;

    System.out.printf(dets+" Average Difference in U %f V %f W %f \n",adifU,adifV,adifW);

    TGCanvas c = new TGCanvas();
    c.setTitle(dets);
    c.draw(hDifU).draw(hDifV,"same").draw(hDifW,"same");
    c.region().showLegend(0.05, 0.95);

    /*TGCanvas c2U = new TGCanvas();
    c2U.setTitle(dets+" U View Position");
    c2U.draw(hDifTrueU);

    TGCanvas c2V = new TGCanvas();
    c2V.setTitle(dets+" V View Position");
    c2V.draw(hDifTrueV);

    TGCanvas c2W = new TGCanvas();
    c2W.setTitle(dets+" W View Position");
    c2W.draw(hDifTrueW);*/

    

}

public void plotCFVarDep(DataList dle,DataList dlv,  EJMLModel model, String dets, int cutVar, 
    String varName, String varUnits,double low, double high,double step) {

        GraphErrors gU = new GraphErrors();
        gU.attr().setMarkerColor(2);
        gU.attr().setMarkerSize(10);
        gU.attr().setTitle("U View");
        gU.attr().setTitleX(varName+" "+varUnits);
        gU.attr().setTitleY(dets+" Position Difference [strips]");

        GraphErrors gV = new GraphErrors();
        gV.attr().setMarkerColor(5);
        gV.attr().setMarkerSize(10);
        gV.attr().setTitle("V View");
        gV.attr().setTitleX(varName+" "+varUnits);
        gV.attr().setTitleY(dets+" Position Difference [strips]");

        GraphErrors gW = new GraphErrors();
        gW.attr().setMarkerColor(3);
        gW.attr().setMarkerSize(10);
        gW.attr().setTitle("W View");
        gW.attr().setTitleX(varName+" "+varUnits);
        gW.attr().setTitleY(dets+" Position Difference [strips]");


        for (double q2=low;q2<high;q2+=step){
            float[] metrics=calcMSE_binned(dle,dlv,model,dets,cutVar,q2,q2+step);
            gU.addPoint(q2+step/2, metrics[0], 0, 0);
            gV.addPoint(q2+step/2, metrics[1], 0, 0);
            gW.addPoint(q2+step/2, metrics[2], 0, 0);
        } // Increment threshold on response

        

        TGCanvas c = new TGCanvas();
        c.setTitle(dets+" Position Difference");
        c.draw(gU).draw(gW,"same").draw(gV,"same");
        c.region().axisLimitsY(-5,5);
        c.region().showLegend(0.6, 0.25);

}

public float[] calcMSE_binned(DataList dle,DataList dlv,  EJMLModel model, String dets,int cutVar,double low,double high){
    int nEvents = dle.getList().size(), nPreds=0;
    float[] output_l = new float[9];

    float adifU=0,adifV=0,adifW=0;

    int out_start=0;
    if(dets=="ECIN"){out_start=3;}
    if(dets=="ECOUT"){out_start=6;}

    for(int i=0;i<nEvents;i+=1) {
        float[] vars = dlv.getList().get(i).floatFirst();
        float[] input = dle.getList().get(i).floatFirst();
        float[] desired_l = dle.getList().get(i).floatSecond();
        model.getOutput(input, output_l);
        for(int j=0;j<9;j++){
            desired_l[j]=desired_l[j]*500;
            output_l[j]=output_l[j]*500;
        }
        int[] desired=convLtoStrip(desired_l);
        int[] output=convLtoStrip(output_l);

        if (vars[cutVar] > low && vars[cutVar]<high) {
            adifU+=(desired[out_start+0]-output[out_start+0]);
            adifV+=(desired[out_start+1]-output[out_start+1]);
            adifW+=(desired[out_start+2]-output[out_start+2]);
            nPreds++;
        }
    }
    adifU=adifU/nPreds;
    adifV=adifV/nPreds;
    adifW=adifW/nPreds;

    return new float[]{adifU,adifV,adifW};

}

int[] convLtoStrip(float[] Ls){
    int[] strips = new int[]{-2,-2,-2,-2,-2,-2,-2,-2,-2};
    DataList conv = DataList.fromCSV("LtoStrip_convTable.csv",
            DataList.range(0,69), DataList.range(0,1));
    int[] nStrips = new int[]{68,63,63,37,37,37,37,37,37};

    for(int i=0;i<9;i++){

        //System.out.printf("det %d \n",i);
        float[] conv_det_view = conv.getList().get(i).floatFirst();
        //System.out.println("got conv");
        for (int str=0;str<(nStrips[i]);str++){
            //in PCAL V, W take distance from other end
            if(i==1 || i==2){
                if(Ls[i]<conv_det_view[str] && Ls[i]>=conv_det_view[str+1]){
                    strips[i]=str+1; //add 1 as the upper limit applies to strip+1 due to inverted order
                }

            } else {
                if(Ls[i]>=conv_det_view[str] && Ls[i]<conv_det_view[str+1]){
                    strips[i]=str;
                }

            }
            
            //System.out.println("converted");
        }
    }

    return strips;
}
