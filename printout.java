HipoReader r = new HipoReader("output_1.h5"); //output_test.h5


CompositeNode tr = new CompositeNode(32100,2,"i",1200);
CompositeNode pt = new CompositeNode(32100,3,"i",1200);
CompositeNode cf = new CompositeNode(32100,11,"i",1200);

Event event = new Event();
int counter = 0;
Bank recpart=r.getBank("REC::Particle");
Bank rectrack=r.getBank("REC::Track");
while(r.hasNext() && counter<100){
    
    r.next(event);
    event.read(tr,32100,2);
    event.read(pt,32100,3);
    event.read(cf,32100,11);
    event.read(recpart);
    event.read(rectrack);
    int nrows = tr.getRows();
    if(nrows>0){
      counter++;
      System.out.println("--- next event\n");
      for(int i = 0; i < nrows; i++){
	      System.out.printf("sector = %2d, charge = %3d, pid = %3d segments [ ", tr.getInt(1,i),tr.getInt(2,i),pt.getInt(3,i));
	      for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
	      System.out.printf("] p %9.5f  %9.5f %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i));
        //System.out.printf("ECAL Strips [ ");
        //for(int k = 1; k < 10; k++) System.out.printf("%9.1f ",cf.getDouble(k,i));
        //System.out.printf("] FTOF Position %9.1f \n\n", cf.getDouble(10,i));
      }
      //recpart.show();
      //rectrack.show();
    }
    

}

/exit
