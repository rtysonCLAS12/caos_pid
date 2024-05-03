HipoReader r = new HipoReader("output_noRecEl_InstaEl_inbending.h5"); //_noRecEl_InstaEl_

System.out.println("read");


CompositeNode tr = new CompositeNode(32100,2,"i",1200);
CompositeNode pt = new CompositeNode(32100,3,"i",1200);
CompositeNode cf = new CompositeNode(32100,11,"i",1200);
CompositeNode pt_test = new CompositeNode(32100,22,"i",1200);
CompositeNode cf_test = new CompositeNode(32100,21,"i",1200);

Event event = new Event();
int counter = 0;
while(r.hasNext()){
    counter++;
    r.next(event);
    event.read(tr,32100,2);
    event.read(pt,32100,3);
    event.read(cf,32100,11);
    event.read(pt_test,32100,22);
    event.read(cf_test,32100,21);

    int nrows = tr.getRows();
    System.out.println("\n--- next event");
    
    int shown=0;
    for(int i = 0; i < nrows; i++){
      if(tr.getInt(1,i)==1){
        if(shown==0){
          System.out.printf("--- InstaRec event %d\n\n",pt.getInt(24,0));
        }
        double Px=pt.getDouble(4,i);
        double Py=pt.getDouble(5,i);
        double Pz=pt.getDouble(6,i);
        double P=Math.sqrt(Px*Px+Py*Py+Pz*Pz);
        double Theta = (float) Math.acos((float)Pz / P)*(180.0 / Math.PI);// Math.atan2(Math.sqrt(Px*Px+Py*Py),Pz);
        double Phi = (float) Math.atan2(Py, Px)*(180.0 / Math.PI);
        System.out.printf("sector = %2d, charge = %3d, pid = %3d, resp= %3.3f segments [ ", tr.getInt(1,i),pt.getInt(1,i),pt.getInt(3,i),pt.getDouble(10,i));
	      for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
	      System.out.printf("] p %9.5f  %9.5f %9.5f, theta %9.5f, phi %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i),Theta,Phi);
        System.out.printf("ECAL Strips [ ");
        for(int k = 1; k < 10; k++) System.out.printf("%9.1f ",cf.getDouble(k,i));
        System.out.printf("] FTOF Position %9.1f\npcal ADC = %9.5f, ecin ADC = %9.5f, ecout ADC = %9.5f, sum HTCC ADC %9.5f \n", cf.getDouble(10,i),pt.getDouble(12,i),pt.getDouble(13,i),pt.getDouble(14,i),pt.getDouble(15,i));
        System.out.printf("HTCC PMTS Left: ");
        for(int k = 16; k < 20; k++) System.out.printf(" %9.5f ",pt.getDouble(k,i));
        System.out.printf("\nHTCC PMTS Right: ");
        for(int k = 20; k < 24; k++) System.out.printf(" %9.5f ",pt.getDouble(k,i));
        System.out.println("\n");
        shown++;
      }
	    
    }
    
    int nrowstest=pt_test.getRows();
    int shown_r=0;
    for(int i = 0; i < nrowstest; i++){
      if(pt_test.getInt(2,i)==1){
        if(shown_r==0){
          System.out.printf("--- REC AI event %d\n\n",pt_test.getInt(17,0));
        }
        double Px=pt_test.getDouble(4,i);
        double Py=pt_test.getDouble(5,i);
        double Pz=pt_test.getDouble(6,i);
        double P=Math.sqrt(Px*Px+Py*Py+Pz*Pz);
        double Theta = (float) Math.acos((float)Pz / P)*(180.0 / Math.PI);// Math.atan2(Math.sqrt(Px*Px+Py*Py),Pz);
        double Phi = (float) Math.atan2(Py, Px)*(180.0 / Math.PI);
	      System.out.printf("sector = %2d, charge = %3d, pid = %3d, status = %3d, track chi^2 = %9.5f ", pt_test.getInt(2,i),pt_test.getInt(1,i),pt_test.getInt(3,i),pt_test.getInt(11,i),pt_test.getDouble(12,i));
	      System.out.printf("p %9.5f  %9.5f %9.5f, theta %9.5f, phi %9.5f\n",pt_test.getDouble(4,i),pt_test.getDouble(5,i),pt_test.getDouble(6,i),Theta,Phi);
        System.out.printf("ECAL Strips [ ");
        for(int k = 1; k < 10; k++) System.out.printf("%9.1f ",cf_test.getDouble(k,i));
        System.out.printf("] FTOF Position %9.1f \n", cf_test.getDouble(10,i));
        System.out.printf("pcal e = %9.5f, ecin e = %9.5f, ecout e = %9.5f, nb nphe = %9.5f \n\n", pt_test.getDouble(13,i), pt_test.getDouble(14,i), pt_test.getDouble(15,i), pt_test.getDouble(16,i));
        shown_r++;
      }
    }
    

}

/exit
