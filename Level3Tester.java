
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
import j4np.neural.classifier.NeuralClassifierModel;

/**
 *
 * @author tyson
 */
public class Level3Tester {
  NeuralClassifierModel trackfinder = new NeuralClassifierModel();
  EJMLModel cf;
  DataList LtoSconv;
  EJMLModel pider;

  List<float[]> negTracks = new ArrayList<>();
  List<float[]> posTracks = new ArrayList<>();
  List<float[]> all_permutations = new ArrayList<>();

  public Level3Tester() {

  }

  public void load_trackfinder(String path) {
    trackfinder.loadFromFile(path, 12);
  }

  public void load_cf(String path) {
    cf = new EJMLModel(path, ModelType.TANH_LINEAR);
  }

  public void load_LtoStripConv(String path) {
    LtoSconv = DataList.fromCSV("LtoStrip_convTable.csv",
        DataList.range(0, 69), DataList.range(0, 1));
  }

  public void load_pider(String path) {
    pider = new EJMLModel(path, ModelType.SOFTMAX);
  }

  public Boolean arrayContains(float[] arr, int c) {
    for (double value : arr) {
      if (value == c) {
        return true;
      }
    }
    return false;
  }

  public void find_tracks_fromClusters(Bank ClusterBank, List<float[]> negTracks, List<float[]> posTracks) {

    all_permutations.clear();

    for (int sect = 1; sect < 7; sect++) {
      List<Float> SL1 = new ArrayList<>();
      List<Float> SL2 = new ArrayList<>();
      List<Float> SL3 = new ArrayList<>();
      List<Float> SL4 = new ArrayList<>();
      List<Float> SL5 = new ArrayList<>();
      List<Float> SL6 = new ArrayList<>();

      for (int i = 0; i < ClusterBank.getRows(); i++) {
        int sl = ClusterBank.getInt("superlayer", i);
        int s = ClusterBank.getInt("sector", i);
        float avgw = ClusterBank.getFloat("avgWire", i) / 112;
        // System.out.println("Adding sls");
        if (s == sect) {
          if (sl == 1) {
            SL1.add(avgw);
          } else if (sl == 2) {
            SL2.add(avgw);
          } else if (sl == 3) {
            SL3.add(avgw);
          } else if (sl == 4) {
            SL4.add(avgw);
          } else if (sl == 5) {
            SL5.add(avgw);
          } else if (sl == 6) {
            SL6.add(avgw);
          }
        }
      }

      // Combine all SLs into a List of Lists
      List<List<Float>> SLs = Arrays.asList(
          SL1, SL2, SL3, SL4, SL5, SL6);

      // Create a list to store permutations
      List<float[]> permutations = new ArrayList<>();

      // System.out.println("gen perms");

      // Call the generatePermutations method to populate the permutations list
      generatePermutations(SLs, 0, new float[6], permutations, 50.0f / 112.0f);

      // for debugging
      for (float[] p : permutations) {
        all_permutations.add(p);
      }

      // System.out.println("Adding perms");
      // Print the permutations
      for (float[] permutation : permutations) {
        float[] tr_pred = new float[3];
        trackfinder.getModel().getOutput(permutation, tr_pred);
        if (tr_pred[2] > 0.1) {
          negTracks.add(permutation);
        } else if (tr_pred[1] > 0.1) {
          posTracks.add(permutation);
        }

      }

    }

  }

  // Recursive method to generate permutations with constraint
  private static void generatePermutations(List<List<Float>> vectors, int depth, float[] current,
      List<float[]> permutations, float maxDifference) {
    if (depth == vectors.size()) {
      // Base case: we've filled in all slots, add this permutation to the list
      permutations.add(Arrays.copyOf(current, current.length));
      return;
    }

    // Iterate over the elements in the current vector
    for (Float element : vectors.get(depth)) {
      // Check the constraint with the previous element
      if (depth > 0 && Math.abs(element - current[depth - 1]) > maxDifference) {
        continue; // Skip this element if the constraint is violated
      }

      // Place the current element in the current slot
      current[depth] = element;

      // Recursively generate permutations for the next depth
      generatePermutations(vectors, depth + 1, current, permutations, maxDifference);
    }
  }

  int[] convLtoStrip(float[] Ls) {
    int[] strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
    int[] nStrips = new int[] { 68, 63, 63, 37, 37, 37, 37, 37, 37 };
    for (int i = 0; i < 9; i++) {
      Ls[i] = Ls[i] * 500;
      float[] conv_det_view = LtoSconv.getList().get(i).floatFirst();
      for (int str = 0; str < (nStrips[i]); str++) {
        // in PCAL V, W take distance from other end
        // System.out.printf("wtf %d %d \n",i,str);
        if (i == 1 || i == 2) {
          if (Ls[i] < conv_det_view[str] && Ls[i] >= conv_det_view[str + 1]) {
            strips[i] = str + 1; // add 1 as the upper limit applies to strip+1 due to inverted order
          }
        } else {
          if (Ls[i] >= conv_det_view[str] && Ls[i] < conv_det_view[str + 1]) {
            strips[i] = str;
          }
        }
      }
    }
    return strips;
  }

  public int nPart_pSect(List<Level3Candidate> Candidates, int sect) {
    int nPart_pSect = 0;
    for (Level3Candidate part : Candidates) {
      if (part.Cal_Sector == sect) {
        nPart_pSect++;
      }
    }
    return nPart_pSect;
  }

  public Boolean arrayNotEmpty(float[] arr) {
    for (double value : arr) {
      if (value != 0) {
        return true;
      }
    }
    return false;
  }

  public int nTrack_pSect(Bank TrackBank, int sect) {
    int nTrack_pSect = 0;
    for (int k = 0; k < TrackBank.getRows(); k++) {
      int pindex = TrackBank.getInt("pindex", k);
      int sectorTrk = TrackBank.getInt("sector", k);
      if (sectorTrk == sect) {
        nTrack_pSect++;
      }
    }
    return nTrack_pSect;
  }

  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc",
  // "HTCC::adc",
  // "REC::Particle", "REC::Calorimeter",
  // "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters");
  public List<Level3Candidate> getCandidates_ownc(Bank[] banks) {

    negTracks.clear();
    posTracks.clear();

    find_tracks_fromClusters(banks[9], negTracks, posTracks);

    // Instant start_time = Instant.now();

    List<Level3Candidate> ps = new ArrayList<Level3Candidate>();
    for (float[] negtrack : negTracks) {
      // find and initialise Candidates
      Level3Candidate part = new Level3Candidate();
      // part.setShow(true);

      part.find_RECParticle_fromtrack(banks[5], banks[8], banks[10], banks[9], negtrack, 0.05);
      part.read_Cal_Bank(banks[6]);
      part.find_HTCC_ADCs(banks[4]);
      part.set_pid_label();

      float[] cf_pred = new float[9];
      cf.getOutput(part.track_clusters, cf_pred);
      int[] cf_strips = convLtoStrip(cf_pred);
      part.read_cal_bank_from_cf_pred(cf_pred, cf_strips, banks[1]);
      float[] pid_pred = new float[2];
      pider.getOutput(part.get_vars_forpid(), pid_pred);
      part.setPidResp(pid_pred[0]);
      // part.print();
      // part.printbanks(tr,pt,banks[5],banks[6],banks[4]);
      ps.add(part);
    }

    // Instant finish_time = Instant.now();
    // timeElapsed += Duration.between(start_time, finish_time).toMillis();

    return ps;
  }

  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc",
  // "HTCC::adc",
  // "REC::Particle", "REC::Calorimeter",
  // "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters");
  public List<Level3Candidate> getCandidates(Bank[] banks, CompositeNode tr, CompositeNode pt) {

    // Instant start_time = Instant.now();

    List<Level3Candidate> ps = new ArrayList<Level3Candidate>();
    for (int i = 0; i < tr.getRows(); i++) {
      // find and initialise Candidates
      Level3Candidate part = new Level3Candidate();
      // part.setShow(true);
      part.find_AI_tracks(tr, pt, i);
      part.find_ClosestRECParticle_fromPredP(banks[5], 999, 999, 999);
      part.read_Cal_Bank(banks[6]);
      part.find_HTCC_ADCs(banks[4]);
      part.set_pid_label();

      float[] cf_pred = new float[9];
      cf.getOutput(part.track_clusters, cf_pred);
      int[] cf_strips = convLtoStrip(cf_pred);
      part.read_cal_bank_from_cf_pred(cf_pred, cf_strips, banks[1]);
      float[] pid_pred = new float[2];
      pider.getOutput(part.get_vars_forpid(), pid_pred);
      part.setPidResp(pid_pred[0]);
      // part.print();
      // part.printbanks(tr,pt,banks[5],banks[6],banks[4]);

      ps.add(part);
    }

    // Instant finish_time = Instant.now();
    // timeElapsed += Duration.between(start_time, finish_time).toMillis();

    return ps;
  }

  public int hasMatchedPID(List<Level3Candidate> Candidates, int pid) {
    int c = 0;
    for (Level3Candidate part : Candidates) {
      int PID = part.PID;
      if (PID == pid) {
        c++;
      }
    }
    return c;
  }

  public int sectHasMatchedPID(List<Level3Candidate> Candidates, int pid, int sect) {
    int c = 0;
    for (Level3Candidate part : Candidates) {
      int PID = part.PID;
      if (PID == pid && part.Sector == sect) {
        c++;
      }
    }
    return c;
  }

  public static int[] convertL1Trigger(long bits) {
    int[] trigger = new int[32];

    // System.out.printf("%X - %X\n", bits,bits&0xF);
    for (int i = 0; i < trigger.length; i++) {
      trigger[i] = 0;
      if (((bits >> i) & (1L)) != 0L)
        trigger[i] = 1;
      // System.out.println(Arrays.toString(trigger));
    }
    return trigger;
  }

  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc",
  // "HTCC::adc",
  // "REC::Particle", "REC::Calorimeter",
  // "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters");
  public int eventHasRECPID(Bank[] banks, int PID, ArrayList<Level3Particle> parts, CompositeNode tr,
      CompositeNode pt) {
    int c = 0;
    for (int i = 0; i < banks[5].getRows(); i++) {
      Level3Particle part = new Level3Particle();
      part.read_Particle_Bank(i, banks[5]);
      part.find_sector_track(banks[8]);
      part.find_track_clusterIDs(banks[10]);
      part.find_track_clusters(banks[9]);
      part.find_HTCC_ADCs(banks[4]);
      part.read_Cal_Bank(banks[6]);
      part.find_sector_cal(banks[6]);

      // System.out.printf("part cal sector %d pcal energy %f
      // \n",part.Cal_Sector,part.PCAL_energy);

      // part.print();

      if (part.Sector != 0 && PID == part.PID && part.Track_nSL == 6) {
        if (PID == 11) {
          if (arrayNotEmpty(part.HTCC_adcs) && part.Vz < 12 && part.Vz > -13) {
            c++;
          }
        } else {
          c++;
        }

      }
    }
    return c;
  }

  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc",
  // "HTCC::adc",
  // "REC::Particle", "REC::Calorimeter",
  // "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters");
  public int sectorHasRECPID(Bank[] banks, int PID, ArrayList<Level3Particle> parts, CompositeNode tr, CompositeNode pt,
      int sect) {
    int c = 0;
    for (int i = 0; i < banks[5].getRows(); i++) {
      Level3Particle part = new Level3Particle();
      part.read_Particle_Bank(i, banks[5]);
      part.find_sector_track(banks[8]);
      part.find_track_clusterIDs(banks[10]);
      part.find_track_clusters(banks[9]);
      part.find_HTCC_ADCs(banks[4]);
      part.read_Cal_Bank(banks[6]);
      part.find_sector_cal(banks[6]);

      // System.out.printf("part cal sector %d pcal energy %f
      // \n",part.Cal_Sector,part.PCAL_energy);

      // part.print();
      // && part.Track_chi2<350 && part.Sector==part.Cal_Sector
      if (part.Sector == sect && PID == part.PID && part.Track_nSL == 6) {
        if (PID == 11) {
          if (arrayNotEmpty(part.HTCC_adcs) && part.Vz < 12 && part.Vz > -13) {
            c++;
          }
        } else {
          c++;
        }

      }
    }
    return c;
  }

  public int sectorHas5SLRECPID(Bank[] banks, int PID, ArrayList<Level3Particle> parts, CompositeNode tr,
      CompositeNode pt, int sect) {
    int c = 0;
    for (int i = 0; i < banks[5].getRows(); i++) {
      Level3Particle part = new Level3Particle();
      part.read_Particle_Bank(i, banks[5]);
      part.find_sector_track(banks[8]);
      part.find_track_clusterIDs(banks[10]);
      part.find_track_clusters(banks[9]);
      part.find_HTCC_ADCs(banks[4]);
      part.read_Cal_Bank(banks[6]);
      part.find_sector_cal(banks[6]);

      // System.out.printf("part cal sector %d pcal energy %f
      // \n",part.Cal_Sector,part.PCAL_energy);

      // part.print();
      // && part.Track_chi2<350 && part.Sector==part.Cal_Sector
      if (part.Sector == sect && PID == part.PID && part.Track_nSL == 5) {
        if (PID == 11) {
          if (arrayNotEmpty(part.HTCC_adcs) && part.Vz < 12 && part.Vz > -13) {
            c++;
          }
        } else {
          c++;
        }

      }
    }
    return c;
  }

  public DataList getPred(String file, int nPart) {

    int nFound = 0, rEv = 0, nRECEls = 0, nREC5SLEls = 0,nRec2Tracks=0, nAI2Tracks=0;
    DataList pred = new DataList();

    // String file=dir;

    HipoReader r = new HipoReader(file);

    Event e = new Event();

    // r.getSchemaFactory().show();
    Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "REC::Particle", "REC::Calorimeter", "REC::Cherenkov", "REC::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters");

    CompositeNode tr = new CompositeNode(3210, 2, "i", 1200);
    CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);

    while (r.hasNext() && nFound < nPart) {
      rEv++;

      r.nextEvent(e);
      e.read(banks);

      // System.out.println("DC");
      // banks[0].show();

      /*
       * System.out.println("FTOF");
       * banks[3].show();
       * System.out.println("HTCC");
       * banks[4].show();
       */

      /*
       * System.out.println("REC::Calorimeter");
       * banks[7].show();
       * System.out.println("ECAL::hits");
       * banks[15].show();
       */

      e.read(tr, 32100, 2);
      e.read(pt, 32100, 3);

      List<Level3Candidate> Candidates = getCandidates(banks, tr, pt);

      for (int sect = 1; sect < 7; sect++) {

        ArrayList<Level3Particle> parts = new ArrayList<Level3Particle>();
        // int nRecEl = eventHasRECPID(banks, 11, parts, tr, pt);
        // int nMatchedEl=hasMatchedPID(Candidates, 11);
        int nRecEl = sectorHasRECPID(banks, 11, parts, tr, pt, sect);
        int nRec5SLEl = sectorHas5SLRECPID(banks, 11, parts, tr, pt, sect);
        int nMatchedEl = sectHasMatchedPID(Candidates, 11, sect);

        nRECEls += nRecEl;
        nREC5SLEls += nRec5SLEl;
        int hasEl = 0;
        if (nRecEl > 0) {
          hasEl = 1;
        }

        int has5SLEl = 0;
        if (nRec5SLEl > 0) {
          has5SLEl = 1;
        }

        double bestp = 0, besttheta = 0, bestphi = -200, bestpid = 0, bestsector = 0;
        double highp = 0, hightheta = 0, highphi = -200, highpid = 0, highsector = 0, maxP = 0;
        int nTracks=0, nAITracks=0;
        for (int i = 0; i < banks[5].getRows(); i++) {
          Level3Particle part = new Level3Particle();
          part.read_Particle_Bank(i, banks[5]);
          part.find_sector_track(banks[8]);
          part.find_HTCC_ADCs(banks[4]);
          // && part.PID!=0 && part.Track_chi2<350 && part.PID!=0&& Math.abs(part.Vz)<20
          if (part.Sector == sect && part.Track_nSL == 6 && part.Charge!=0 && part.PID!=Math.abs(11) && part.Track_chi2<350  && part.PID!=0&& Math.abs(part.Vz)<20) {
            nTracks++;
          }
        }

        /*for (int i = 0; i < banks[8].getRows(); i++) {
          double chi2 = banks[8].getFloat("chi2",i);
          int sectorTrk = banks[8].getInt("sector", i);
          if( (banks[8].getInt("status", i) & 0b101010101010)==0){
            if( chi2<350 &&sectorTrk==sect){ 
              nTracks++;
            }
          }
        }*/

        /*for(int i=0;i<tr.getRows();i++){
          if(tr.getInt(1,i)==sect){
            nAITracks++;  
          }
        }*/

        float maxPred = 0;
        for (Level3Candidate p : Candidates) {
          if ( p.Sector == sect) {
            if(p.unmatched == false && p.Charge == -1){
              if (p.pid_resp > maxPred) {
                maxPred = p.pid_resp;
                bestp = p.Pred_P;
                bestsector = p.Sector;
                bestphi = p.Pred_Phi * (180.0 / Math.PI);
                besttheta = p.Pred_Theta * (180.0 / Math.PI);
                bestpid = p.PID;
              }
            }

            if( p.nCalHits()>-1 ){ //&& p.unmatched == false
              if(p.Charge==-1 && p.pid_resp<0.1){
                nAITracks++;
                if (p.P > maxP) {
                  highp = p.Pred_P;
                  highsector = p.Sector;
                  highphi = p.Pred_Phi * (180.0 / Math.PI);
                  hightheta = p.Pred_Theta * (180.0 / Math.PI);
                  highpid = p.PID;
                  maxP = p.Pred_P;
                }
              } else{
                nAITracks++;
                if (p.P > maxP) {
                  highp = p.Pred_P;
                  highsector = p.Sector;
                  highphi = p.Pred_Phi * (180.0 / Math.PI);
                  hightheta = p.Pred_Theta * (180.0 / Math.PI);
                  highpid = p.PID;
                  maxP = p.Pred_P;
                }
              }
            }
          }
        }

        long bits = banks[2].getLong("trigger", 0);
        int[] L1trigger = convertL1Trigger(bits);

        // don't want to complicate metrics by change increase in efficiency of AI vs
        // REC
        if (nMatchedEl >= nRecEl) {
          if (hasEl == 0 && has5SLEl == 1) {
            // L1 trigger might fire for 5 SL electron
            // but we're only calculating metrics for 6SL electron
            // don't want to unfairly decrease L1 trigger purity
          } else {
            int pidout = 0;
            if (hasEl == 1) {
              pidout = 11;
            }

            int twoTracks=0,twoAITracks=0;
            if (nTracks>1){
              twoTracks=1;
              nRec2Tracks++;
            }
            if (nAITracks>1){
              twoAITracks=1;
              nAI2Tracks++;
            }

            if(nAITracks>1 && nTracks<2){
              /*banks[5].show();
              banks[8].show();
              banks[9].show();
              banks[10].show();
              for(int i=0;i<tr.getRows();i++){
                System.out.printf("sector = %2d, charge = %3d segments [ ", tr.getInt(1,i),tr.getInt(2,i));
                for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
                System.out.printf("] p %9.5f  %9.5f %9.5f %9.5f \n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i),Math.sqrt(pt.getDouble(4,i)*pt.getDouble(4,i)+pt.getDouble(5,i)*pt.getDouble(5,i)+pt.getDouble(6,i)*pt.getDouble(6,i))); 
              }*/
            }

            DataEntry de = new DataEntry(
                new double[] { pidout, bestp, besttheta, bestphi, sect, L1trigger[0 + sect], L1trigger[14 + sect], twoTracks,highp}, // L1trigger[0+sect]
                new double[] { maxPred, twoAITracks });
            pred.add(de);
            nFound++;

          }

        }

      }

    }

    System.out.printf("number of REC e- %d and 5SL e- %d \n", nRECEls, nREC5SLEls);
    System.out.printf("number of sector with 2 REC tracks %d, 2 AI tracks %d \n", nRec2Tracks, nAI2Tracks);
    System.out.printf("Read %d events\n\n", rEv);

    return pred;

  }

  public void PlotResponse(DataList pred, int elClass, int desiredPID, String part) {
    int NEvents = pred.getList().size();

    H1F hRespPos = new H1F(part + " in Sector", 100, 0, 1);
    hRespPos.attr().setLineColor(2);
    hRespPos.attr().setFillColor(2);
    hRespPos.attr().setLineWidth(3);
    hRespPos.attr().setTitleX("Response");
    H1F hRespNeg = new H1F("No " + part + " in Sector", 100, 0, 1);
    hRespNeg.attr().setLineColor(5);
    hRespNeg.attr().setLineWidth(3);
    hRespNeg.attr().setTitleX("Response");
    // Sort predictions into those made on the positive/or negative samples
    for (int i = 0; i < NEvents; i += 1) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] resp = pred.getList().get(i).floatSecond();

      if (vars[0] == desiredPID) {
        hRespPos.fill(resp[elClass]);
      } else {
        hRespNeg.fill(resp[elClass]);
      }
    }

    TGCanvas c = new TGCanvas();

    c.setTitle("Response");
    c.draw(hRespPos).draw(hRespNeg, "same");
    c.region().showLegend(0.05, 0.95);

  }// End of PlotResponse

  public void PlotVar(DataList pred, int varIndex, int desiredPID, String part, String varName, String varUnits,
      double low, double high, int nBins, double respTh, int elClass) {
    int NEvents = pred.getList().size();

    int noPlot = 0;
    if (varIndex == 3) {
      noPlot = -5;
    }

    H1F hRespPos = new H1F(part + " in Sector", nBins, low, high);
    hRespPos.attr().setLineColor(2);
    hRespPos.attr().setFillColor(2);
    hRespPos.attr().setLineWidth(3);
    hRespPos.attr().setTitleX(varName + " " + varUnits);
    H1F hRespNeg = new H1F("No " + part + " in Sector", nBins, low, high);
    hRespNeg.attr().setLineColor(5);
    hRespNeg.attr().setLineWidth(3);
    hRespNeg.attr().setTitleX(varName + " " + varUnits);
    // Sort predictions into those made on the positive/or negative samples
    for (int i = 0; i < NEvents; i += 1) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] resp = pred.getList().get(i).floatSecond();

      if (resp[elClass] > respTh) {// && vars[varIndex]!=noPlot
        if (vars[0] == desiredPID) {
          hRespPos.fill(vars[varIndex]);
        } else {
          hRespNeg.fill(vars[varIndex]);
        }
      }
    }

    TGCanvas c = new TGCanvas();

    c.setTitle(varName);
    c.draw(hRespPos).draw(hRespNeg, "same");
    c.region().showLegend(0.05, 0.95);

  }// End of PlotResponse

  public StatNumber[] getMetrics(DataList pred, int elClass, int desiredPID, double thresh) {

    int nEvents = pred.getList().size();

    double TP = 0, FP = 0, FN = 0;
    for (int i = 0; i < nEvents; i++) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] resp = pred.getList().get(i).floatSecond();
      if (vars[0] == desiredPID) {
        if (resp[elClass] > thresh) {
          TP++;
        } else {
          FN++;
        }
      } else {
        if (resp[elClass] > thresh) {
          FP++;
        }
      } // Check true label
    }

    /*
     * System.out.printf("Theres %d electrons in sample\n", nEls);
     * System.out.printf("L1 trigger fired %d times in sample\n", nTrig);
     */
    return calcMetrics(TP, FP, FN);
  }

  public StatNumber[] calcMetrics(double TP, double FP, double FN) {
    StatNumber Pur = new StatNumber(TP, Math.sqrt(TP));
    StatNumber Eff = new StatNumber(TP, Math.sqrt(TP));
    StatNumber FPs = new StatNumber(FP, Math.sqrt(FP));
    StatNumber FNs = new StatNumber(FN, Math.sqrt(FN));
    StatNumber denomPur = new StatNumber(TP, Math.sqrt(TP));
    StatNumber denomEff = new StatNumber(TP, Math.sqrt(TP));
    denomPur.add(FPs);
    denomEff.add(FNs);
    Pur.divide(denomPur);
    Eff.divide(denomEff);
    StatNumber[] mets = new StatNumber[2];
    mets[0] = Pur;
    mets[1] = Eff;
    return mets;
  }

  // Labels col 0 is 1 if there's an e-, 0 otherwise
  public StatNumber[] getMetsForBin(DataList pred, int elClass, int targetVar, int L1ind, int desiredPID, double thresh, int cutVar,
      double low, double high) {

    int nEvents = pred.getList().size();

    //L1ind is -1 we're ignoreing L1 trigger, so can just take same number as L3
    if(L1ind==-1){L1ind=targetVar;}

    double TP = 0, FN = 0, FP = 0;
    double TPL1 = 0, FNL1 = 0, FPL1 = 0;
    for (int i = 0; i < nEvents; i++) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] resp = pred.getList().get(i).floatSecond();
      if (vars[cutVar] > low && vars[cutVar] < high) {
        if (vars[targetVar] == desiredPID) {
          if (resp[elClass] > thresh) {
            TP++;
          } else {
            FN++;
          }
          if (vars[L1ind] == 1) {
            TPL1++;
          } else {
            FNL1++;
          }
        } else {
          if (resp[elClass] > thresh) {
            FP++;
          }
          if (vars[L1ind] == 1) {
            FPL1++;
          }
        } // Check true label
      }
    }
    StatNumber[] mets = new StatNumber[4];
    StatNumber[] metsL3 = calcMetrics(TP, FP, FN);
    mets[0] = metsL3[0];
    mets[1] = metsL3[1];

    StatNumber[] metsL1 = calcMetrics(TPL1, FPL1, FNL1);
    mets[2] = metsL1[0];
    mets[3] = metsL1[1];

    return mets;
  }

  public void plotVarDep(DataList pred, int elClass, int targetVar, int desiredPID, double thresh, Boolean addPur,
      int L1ind, int cutVar, String varName, String varUnits, double low, double high, double step) {

    String yTitle = "Metrics";
    if (!addPur) {
      yTitle = "Efficiency";
    }

    GraphErrors gEff = new GraphErrors();
    gEff.attr().setMarkerColor(2);
    gEff.attr().setMarkerSize(10);
    gEff.attr().setTitle("Level3 Efficiency");
    gEff.attr().setTitleX(varName + " " + varUnits);
    gEff.attr().setTitleY(yTitle);

    GraphErrors gPur = new GraphErrors();
    gPur.attr().setMarkerColor(5);
    gPur.attr().setMarkerSize(10);
    gPur.attr().setTitle("Level3 Purity");
    gPur.attr().setTitleX(varName + " " + varUnits);
    gPur.attr().setTitleY(yTitle);

    GraphErrors gEffL1 = new GraphErrors();
    gEffL1.attr().setMarkerColor(3);
    gEffL1.attr().setMarkerSize(10);
    gEffL1.attr().setTitle("Level1 Efficiency");
    gEffL1.attr().setTitleX(varName + " " + varUnits);
    gEffL1.attr().setTitleY(yTitle);

    GraphErrors gPurL1 = new GraphErrors();
    gPurL1.attr().setMarkerColor(6);
    gPurL1.attr().setMarkerSize(10);
    gPurL1.attr().setTitle("Level1 Purity");
    gPurL1.attr().setTitleX(varName + " " + varUnits);
    gPurL1.attr().setTitleY(yTitle);

    for (double q2 = low; q2 < high; q2 += step) {
      StatNumber[] metrics = getMetsForBin(pred, elClass,targetVar, L1ind, desiredPID, thresh, cutVar, q2, q2 + step);
      gPur.addPoint(q2 + step / 2, metrics[0].number(), 0, metrics[0].error());
      gEff.addPoint(q2 + step / 2, metrics[1].number(), 0, metrics[1].error());
      gPurL1.addPoint(q2 + step / 2, metrics[2].number(), 0, metrics[2].error());
      gEffL1.addPoint(q2 + step / 2, metrics[3].number(), 0, metrics[3].error());
    } // Increment threshold on response

    TGCanvas c = new TGCanvas();
    c.setTitle("Efficiency vs " + varName);
    c.draw(gEff);
    if (addPur) {
      c.draw(gPur, "same");
    }
    if (L1ind != -1) {
      // don't draw L1 efficiency, this is basically 1
      // any divergence from 1 is due to cuts in DC roads trigger
      // misleading to add to plots .draw(gEffL1,"same")
      c.draw(gPurL1, "same");
      c.region().axisLimitsY(gPurL1.getVectorY().getMin() - 0.1, 1.05);
    } else {
      if(addPur){
        c.region().axisLimitsY(gPur.getVectorY().getMin() - 0.1, 1.05);
      } else {
        c.region().axisLimitsY(gEff.getVectorY().getMin() - 0.1, 1.05);
      }
      
    }

    c.region().showLegend(0.6, 0.25);

  }

  public double findBestThreshold(DataList pred, int elClass, int desiredPID, double effLow) {

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
    double bestPuratEffLow = 0;
    double bestPurErratEffLow = 0;

    // Loop over threshold on the response
    for (double RespTh = 0.01; RespTh < 0.99; RespTh += 0.01) {
      StatNumber[] metrics = getMetrics(pred, elClass, desiredPID, RespTh);
      double Pur = metrics[0].number();
      double Eff = metrics[1].number();
      double PurErr = metrics[0].error();
      double EffErr = metrics[1].error();
      gPur.addPoint(RespTh, Pur, 0, PurErr);
      gEff.addPoint(RespTh, Eff, 0, EffErr);
      if (Eff > effLow) {
        if (Pur > bestPuratEffLow) {
          bestPuratEffLow = Pur;
          bestPurErratEffLow = PurErr;
          bestRespTh = RespTh;
        }
      }
    } // Increment threshold on response

    System.out.format(
        "%n Best Purity at Efficiency above %.3f: %.3f +/- %.4f at a threshold on the response of %.3f %n%n",
        effLow, bestPuratEffLow, bestPurErratEffLow, bestRespTh);

    TGCanvas c = new TGCanvas();
    c.setTitle("Metrics vs Response");
    c.draw(gEff).draw(gPur, "same");
    c.region().showLegend(0.25, 0.25);
    // c.region().axisLimitsY(0.8, 1.01);

    return bestRespTh;
  }

  public static void main(String[] args) {

    // to run, in order
    // /open Level3Particle.java
    // /open Level3Candidate.java
    // /open Level3Tester.java
    // Level3Tester.main(new String[]{});

    String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/recook_caos_pid/run_018326_4_wAIBanks.h5";

    Level3Tester tester = new Level3Tester();
    tester.load_trackfinder("clas12rgd.network");
    tester.load_cf("cf_el.network");
    tester.load_LtoStripConv("LtoStrip_convTable.csv");
    tester.load_pider("pid_elNegBG_fromcfpred.network");// _fromcfpred

    DataList preds = tester.getPred(file, 1000000);

    tester.PlotVar(preds, 1, 11, "e-", "P", "[GeV]", 0, 10, 100, -1, 0);
    tester.PlotVar(preds, 2, 11, "e-", "Theta", "[Deg]", 0, 50, 100, -1, 0);
    tester.PlotVar(preds, 3, 11, "e-", "Phi", "[Deg]", -180, 180, 100, -1, 0);

    tester.PlotResponse(preds, 0, 11, "e-");
    double bestth = tester.findBestThreshold(preds, 0, 11, 0.995);
    // set L1ind to -1 to avoid plotting L1 trigger
    tester.plotVarDep(preds, 0,0, 11, bestth, true, 5, 1, "P", "[GeV]", 1, 9.0, 1.0);
    tester.plotVarDep(preds, 0,0, 11, bestth, true, 5, 2, "Theta", "[Deg]", 5.0, 35.0, 5.);
    tester.plotVarDep(preds, 0,0, 11, bestth, true, 5, 3, "Phi", "[Deg]", -180, 180, 10.);
    tester.plotVarDep(preds, 0,0, 11, bestth, true, 5, 4, "Sector", "", 0.5, 6.5, 1.0);


    tester.plotVarDep(preds, 1,7, 1, 0.5, false, -1, 8, "Highest Track P", "[GeV]", 1, 9.0, 1.0);
    tester.plotVarDep(preds, 1,7, 1, 0.5, false, -1, 4, "Sector", "", 0.5, 6.5, 1.0);

  }
}
