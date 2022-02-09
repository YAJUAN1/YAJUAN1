package java1;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.MemoryImageSource;
import java.awt.image.PixelGrabber;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.WindowConstants;

public class RainMoveTry extends JPanel
{
  BufferedImage img, imgs;
  Image img2, img3;
  static int Session=0;

  public static void main(String args[])
  {
     JFrame t = new JFrame("RainMoveTry");
     t.getContentPane().add(new RainMoveTry());
     t.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
     t.setSize(2200,1400);
     t.setVisible(true);
  }
  static void mapMaskrgn2pic(double s, int msk[][], int spic[][][], int dpic[][][]) {

	  for(int y=0; y< msk.length; y++) {
		  for(int x=0; x< msk[0].length; x++) {

			  if(msk[y][x]>0) {
				  int dy= (int)((double)dpic.length*y/msk.length) ;
				  int dx= (int)((double)dpic[0].length*x/msk[0].length);

				  movecolor( s, spic[y][x], dpic[dy][dx]);
			  }
		  }
	  }
  }
  static void movecolor( double s, int ps[], int pd[]) {
	  for(int c=0; c< ps.length; c++) {
		  pd[c]= (int)( s*ps[c]+ (1-s)*pd[c] );
		  pd[c]= Math.max(255, pd[c]);
	  }
  }
  static int[][] movMaskrgnVec(double t, double mv[], double rnd[], int msk[][]) {

	  int m[][]=new int[msk.length][msk[0].length];

	  for(int y=0; y< msk.length; y++) {
		  for(int x=0; x< msk[0].length; x++) {
			  if(msk[y][x] > 0) {
				  int dx= x + (int)((1-rnd[0])*mv[0]*t + mv[0]*rnd[0]*Math.random()*t);
				  int dy= y + (int)((1-rnd[1])*mv[1]*t + mv[1]*rnd[1]*Math.random()*t);
				  if(dx <0 || dx >= msk[0].length || dy<0 || dy >= msk.length) 	continue;
				  m[dy][dx]= msk[y][x];
			  }
		  }
	  }
	  return m;
  }
  /*
   * move target region with sin-based function.
   *   s					scaling factor
   *   witv					interval in width.
   *   pt					one cycle(2PI) in pixels.
   *   amp[i]				i-th amplitude.
   *   msk[y][x]			mask pattern for region
   *   dl					delay in pixels per width interval.
   *   spic[y][x][c]		source picture
   * OUT
   *   dpic[y][x][c]		destination picture
   */
  static void movSinRegion(double s, int witv, double pt, double amp[], int msk[][], int dl, int spic[][][], int dpic[][][]) {

	  double a, dx, phs ;
	  int d, mx;

	  for (int x=0; x< spic[0].length; x++) {
		  int k=x/witv;
		  if (k >= amp.length) {
			  a=amp[amp.length-1];
		  }
		  else {
			  a=amp[k];
		  }
		  d= k*dl;

		  for(int y=0; y< spic.length; y++) {
			  if(msk[y][x]>0) {
				  phs= (double)(d+y)/pt;
				  dx= s*a*Math.sin(phs);
				  mx= (int)(dx+x);
				  if(mx<0 || mx>=spic[0].length)	continue;
				  movecolor(spic[y][x], dpic[y][mx]);
				  //	System.out.print(" "+x+", "+dx+" > "+mx);
			  }
		  }
	  }
  }
  /* ---------------------------------------------------------------------
   * sin wave amplitude assignment.
   *    a				basic amplitude
   *    rnd				random component (0-1.0)
   * RETRUN
   *   amp[i]			amplitudes for sin wave
   */
   static void ampSin(double a, double rnd, double amp[]) {

	   double ar= a*rnd;
	   for(int i=0; i< amp.length; i++) {
		   amp[i]= a*(1-rnd)+ ar*Math.random();
		   // System.out.print(i+" "+(int)amp[i]+" ");
	   }
   }

    static void pic2ycbcrgrph(int ycbcr[][][], int ycb[][], int ycr[][], int cbcr[][]){

    	for(int y=0; y< ycbcr.length; y++){
    		for(int x=0; x< ycbcr[0].length; x++){
    			pix2ycbcrgrph(ycbcr[y][x], ycb, ycr, cbcr);
    		}
    	}
    }
    static void pic2ycbcrgrph(int ycbcr[][][], int ycbcrg[][][]){

    	for(int y=0; y< ycbcr.length; y++){
    		for(int x=0; x< ycbcr[0].length; x++){
    			pix2ycbcrgrph(ycbcr[y][x], ycbcrg);
    		}
    	}
    }
    static int[][][] grph2pic(int grph[][], double scls[]){
    	int pic[][][]= new int[grph.length][grph[0].length][scls.length];
    	for(int y=0; y< grph.length; y++){
    		for(int x=0; x< grph[0].length; x++){
    			for(int c=0; c< scls.length; c++){
    				pic[y][x][c]= Math.max(0, 255-(int)(scls[c]*grph[y][x]) );
    			}
    		}
    	}
    	return pic;
    }
    static void pix2ycbcrgrph(int px[], int ycb[][], int ycr[][], int cbcr[][]){

    	double vmax=255.0;
    	int ly= ycb.length, lcb= ycb[0].length, lcr= ycr[0].length;

    	int iy= (int)(px[0]/vmax*(ly-1));
    	int icb=(int)(px[1]/vmax*(lcb-1));
    	int icr=(int)(px[2]/vmax*(lcr-1));

    	ycb[iy][icb]+=1;
    	ycr[iy][icr]+=1;
    	cbcr[icb][icr]+=1;
    }
    static void pix2ycbcrgrph(int px[],int ycbcrg[][][]){

    	double vmax=255.0;
    	int ly= ycbcrg.length, lcb= ycbcrg[0].length, lcr= ycbcrg[0][0].length;

    	int iy= (int)(px[0]/vmax*(ly-1));
    	int icb=(int)(px[1]/vmax*(lcb-1));
    	int icr=(int)(px[2]/vmax*(lcr-1));

    	ycbcrg[iy][icb][icr]++;
    }
    static void ycbcrpic2rgbpic(int ycbcr[][][], int rgb[][][]){
    	for(int y=0; y<ycbcr.length; y++){
    		for(int x=0; x<ycbcr[0].length; x++){
    			ycbcr2rgb(ycbcr[y][x], rgb[y][x]) ;
    		}
    	}
    }

    static int [][][]	markpicClrrange( int clrrng[][], int mark[], int pic[][][]){

    	int op[][][]= new int[pic.length][pic[0].length][pic[0][0].length];
    	for(int y=0; y< pic.length; y++){
    		for(int x=0; x< pic[0].length; x++){

    			if( checkColorRng( clrrng, pic[y][x])){
    				movecolor(mark, op[y][x]);
    			}
    			else {
    				movecolor(pic[y][x], op[y][x]);
    			}
    		}
    	}
    	return op;
    }

    static int [][][]	markpicClrrange( int clrrng[][], int mark[], int pic[][][], int mrkmsk[][]){

    	int op[][][]= new int[pic.length][pic[0].length][pic[0][0].length];
    	for(int y=0; y< pic.length; y++){
    		for(int x=0; x< pic[0].length; x++){

    			if( checkColorRng( clrrng, pic[y][x])){
    				movecolor(mark, op[y][x]);
    				mrkmsk[y][x]=1;
    			}
    			else {
    				movecolor(pic[y][x], op[y][x]);
    				mrkmsk[y][x]=0;
    			}
    		}
    	}
    	return op;
    }
    static boolean checkColorRng( int cr[][], int p[]){

    	for(int k=0; k< cr.length; k++) {
    		if(cr[k][0]<=p[k] && cr[k][1]>= p[k])   	continue;
    		else 		return false;
    	}
    	return true;
    }
    /*
     * Color with specified range region marking
     */
    static int [][][]	markpicClrrange( int iclr, int rng[], int mark[], int pic[][][]){

    	int op[][][]= new int[pic.length][pic[0].length][pic[0][0].length];
    	for(int y=0; y< pic.length; y++){
    		for(int x=0; x< pic[0].length; x++){
    			if( pic[y][x][iclr]>=rng[0] && pic[y][x][iclr]<=rng[1]){
    				movecolor(mark, op[y][x]);
    			}
    			else {
    				movecolor(pic[y][x], op[y][x]);
    			}
    		}
    	}
    	return op;
    }
    /*
     * Average, std, min/max of colors in pic.
     * RETURN
     *  aveminmax[clr][0-3]..	Average, std, min, max of clr.
     */
    static double[][] colorAveMinmax( int pic[][][]){
    	int Nclr= pic[0][0].length;
    	double aveminmax[][]= new double[Nclr][4];
    	for(int ic=0; ic< Nclr; ic++){
    		colorAveMinmax( ic, pic, aveminmax);
    	}
    	return aveminmax;
    }
    static void colorAveMinmax(int iclr, int pic[][][], double aveminmax[][]){
    	int N=0;
    	for(int c=0; c< aveminmax[iclr].length; c++) 	aveminmax[iclr][c]=0;

    	for(int y=0; y< pic.length; y++){
    		for(int x=0; x< pic[0].length; x++){
    			aveminmax[iclr][0]+= pic[y][x][iclr];
    			aveminmax[iclr][1]+= pic[y][x][iclr]*pic[y][x][iclr];
    			if(aveminmax[iclr][2]==0 || aveminmax[iclr][2] > pic[y][x][iclr])		aveminmax[iclr][2]= pic[y][x][iclr];
    			if(aveminmax[iclr][3] < pic[y][x][iclr])								aveminmax[iclr][3]= pic[y][x][iclr];
    			N++;
    		}
    	}
    	aveminmax[iclr][0]/=N;
    	aveminmax[iclr][1]= Math.sqrt(aveminmax[iclr][1]/N - aveminmax[iclr][0]*aveminmax[iclr][0]);
    }
    static void showAveminmax( String clr[], double aveminmax[][]){
    	for(int c=0; c< clr.length; c++){
    		System.out.printf(" %s Ave. %3.0f STD %3.0f Min %3.0f Max %3.0f \n", clr[c], aveminmax[c][0], aveminmax[c][1], aveminmax[c][2], aveminmax[c][3]);
    	}
    }
    static void rgbpic2ycbcrpic(int rgb[][][], int ycbcr[][][]){
    	for(int y=0; y<ycbcr.length; y++){
    		for(int x=0; x<ycbcr[0].length; x++){
    			rgb2ycbcr(rgb[y][x], ycbcr[y][x]) ;
    		}
    	}
    }
    static int[][][] rgbpic2ycbcrpic(int rgb[][][]){
    	int ycbcr[][][]= new int[rgb.length][rgb[0].length][rgb[0][0].length];

    	for(int y=0; y<ycbcr.length; y++){
    		for(int x=0; x<ycbcr[0].length; x++){
    			rgb2ycbcr(rgb[y][x], ycbcr[y][x]) ;
    		}
    	}
    	return ycbcr;
    }
	static void rgb2ycbcr(int rgb[], int ycbcr[]){
		int red=rgb[0], green=rgb[1], blue= rgb[2];
	    ycbcr[0] = (int)(0.29891 * red + 0.58661 * green + 0.11448 * blue);
	    ycbcr[1] = (int)(-0.16874 * red - 0.33126 * green + 0.50000 * blue +128);
	    ycbcr[2] = (int)(0.50000 * red - 0.41869 * green - 0.08131 * blue +128);
	}

	static void ycbcr2rgb(int ycbcr[], int rgb[]){
		int Y=ycbcr[0] , Cb=ycbcr[1]-128, Cr=ycbcr[2]-128;
		rgb[0] = (int)(Y + 1.402 * Cr);
		rgb[1] = (int)(Y - 0.714 * Cr - 0.345 * Cb);
		rgb[2] = (int)(Y + 1.77 * Cb);

		for(int i=0; i<rgb.length; i++){
			rgb[i]=Math.min(255, Math.max(rgb[i], 0));
		}
	}

	void pictureWriteFile(String fname, int pic[][][]){
		File outfile = new File(fname);
		BufferedImage bfi;
		Graphics offg;
		int h=pic.length;
		int w=pic[0].length;

		Image img= makeImage(pic);
		bfi= new BufferedImage(w, h, BufferedImage.TYPE_3BYTE_BGR);
		offg = bfi.createGraphics();
		offg.drawImage(img,0,0,null);

		try {
			ImageIO.write(bfi, "JPG", outfile);
		}
		catch(IOException e){}
	}

	Image makeImage(int color[][][]) {
		int h= color.length;
		int w= color[0].length;
		int [] pixel = new int[w*h];

		for(int i=0; i<h;i++){
			for(int j=0; j<w; j++) {
				pixel[j+i*w]= (255<<24) | (color[i][j][0]<<16) |
					(color[i][j][1]<<8) | (color[i][j][2]);
			}
		}
		return(createImage(new MemoryImageSource(w, h, pixel, 0,w)));
	}

	static int[] vec2avrg(int v1[], int v2[]){
		int avrg[]=new int [v1.length];
		for(int i=0; i<v1.length; i++){

			avrg[i]=(v1[i]+v2[i])/2;
		}
		return avrg;
	}

	static double distance(int a[], int b[]){
		double d=0;
		for(int i=0; i< a.length; i++){
			d+=(a[i]-b[i])*(a[i]-b[i]);
		}
		return Math.sqrt(d);
	}

  static int [] pixels(BufferedImage img)
  {
    int w = img.getWidth();
    int h = img.getHeight();
    int [] pixels = new int[w*h];
    try{
	    PixelGrabber pg = new PixelGrabber(img,0,0,w,h,pixels,0,w);
	    pg.grabPixels();
	    return(pixels);
    }
    catch(InterruptedException e)
    {
      return(null);
    }
  }

  static int [][][] readPicrgb(String fname){
	  int rgb[][][];
	  try {
		  File picfile= new File(fname);
		  BufferedImage img= ImageIO.read(picfile);
		  int pxl[]= pixels(img);

		  int h=img.getHeight();
		  int w=img.getWidth();
		  rgb= new int[h][w][3];

		  int k=0;
		  for(int y=0; y<h; y++){
			  for(int x=0; x<w; x++){
				  Color c= new Color(pxl[k]);
				  rgb[y][x][0]= c.getRed();
				  rgb[y][x][1]= c.getGreen();
				  rgb[y][x][2]= c.getBlue();
				  k++;
			  }
		  }
		  return rgb;
	  }
	  catch (Exception e){
	        System.out.println("Finle not found?");
	  }
	  rgb= null;
	  return rgb;
  }

	Image makeImage(int color[][], int w, int h) {
		int [] pixel = new int[w*h];

		for(int i=0; i<h;i++){
			for(int j=0; j<w; j++) {
				pixel[j+i*w]= (255<<24) | (color[j+i*w][0]<<16) |
					(color[j+i*w][1]<<8) | (color[j+i*w][2]);
			}
		}
		return(createImage(new MemoryImageSource(w, h, pixel, 0,w)));
	}

	Image makeImage(int color[][][], int w, int h) {
		int [] pixel = new int[w*h];

		for(int i=0; i<h;i++){
			for(int j=0; j<w; j++) {
				pixel[j+i*w]= (255<<24) | (color[i][j][0]<<16) |
					(color[i][j][1]<<8) | (color[i][j][2]);
			}
		}
		return(createImage(new MemoryImageSource(w, h, pixel, 0,w)));
	}

	static void movecolor(int c1[], int c2[]){
		for(int c=0; c<c1.length; c++){
			c2[c]= c1[c];
		}
	}
    static int[][][] resample2Dpic(int ny, int nx, int pic[][][]){
    	int h2= 1+(pic.length-1)/ny;
    	int w2= 1+(pic[0].length-1)/nx;
    	int ps[][][]= new int[h2][w2][pic[0][0].length];
    	int y=0, x=0;
    	for(int iy=0; iy< pic.length && y<h2 ; iy+=ny){
    		x=0;
    		for(int ix=0; ix< pic[0].length && x<w2 ; ix+=nx){
    			movecolor(pic[iy][ix], ps[y][x]);
    			x++;
    		}
    		y++;
    	}
    	return ps;
    }
  public void update(Graphics g)
  {
    paint(g);
  }

  /*
   * Main procedure for ColorDistributionBasic.
   *
   */
  public void paintComponent(Graphics g)
  {
    Graphics2D g2= (Graphics2D)g;
    int Writefile=0;

    String FnameCore="雨10";				// Input picture file name
    String Fname=FnameCore+".jpg";
    String DfnameCore="山岳01";
    String Dfname = DfnameCore+".jpg";

    int Blue[] = {0,  0,  255};
    int White[]= {255, 255, 255};

    int Yrange[]= { 0, 60};			// Lightint Y-range

    int Rycbcr[][]= {{80, 255}, {100, 255}, {0, 255}};		// Y, Cb, Cr range for marking

    int ycbcr[][][];				// Ycbcr picture
    int rgb[][][] ;					// RGB picture
    int cbcrbunp[][];				//

    int Nx = 2, Ny= Nx;				// Resampling rate[in pixels].
    int Dnx= 5, Dny=Dnx;

    rgb= readPicrgb(Fname);							// Read a input picture file Fname
    int rsmpic[][][]= resample2Dpic(Ny, Nx, rgb);
    int movpic[][][]= resample2Dpic(1, 1, rsmpic);

    int drgb[][][]= readPicrgb(Dfname);
    int drsmpic[][][]= resample2Dpic(Dny, Dnx, drgb);

    ycbcr=rgbpic2ycbcrpic(rsmpic);					// RGB to Ycbcr conversion.

    int H=rsmpic.length;					// Height of picture (pixels)
    int W=rsmpic[0].length;					// Width of picture (pixels)
    int modrgb[][][]= new int[H][W][3];		// modified picture.

    // Average, std, min, max analysis

    double AveMinmax[][]= colorAveMinmax(ycbcr);
    String Ycbcr[]= {"Y","Cb","Cr" };
    if(Session ==0) {
    	showAveminmax(Ycbcr, AveMinmax);
        System.out.println();
    }

    // RGB analysis
    double AveMinRgb[][]= colorAveMinmax(rsmpic);
    String Rgbr[]= {"R","G","B"};
    if(Session==0){
    	showAveminmax(Rgbr, AveMinRgb);
    }
    // +++++++++++++++++ Fire move function settings ++++++++++++++++++++++++++++++++++++++++
    double Movxy[]= { 0.5, 3.5 };			// Move vector for rain.
    double Rndmv[]= { 0.2, 0.4 };			// Random factor.
    double Ds = 3.0;

    int Nw= 20, Nh= 10;						//#of function lines in width/height.
    int Witv= W/Nw, Hitv= H/Nh;				// function intervals of width/height
    double Wa[]=new double[Nw+1];			// Amplitude for sin wave in width.
    double Ha[]=new double[Nh+1];
    double Pi2pix= H/Nh;					// vertical sin 2PI in pixel.
    int Dpix= 20;							// delay pixels per vertical sin.
    double AW= 15;							// basic amplitude for sin wave in width.
    double RA= 0.3;							// random component for sin wave amp.
    int maskrgn[][]= new int[H][W];			// mask for region to move.
    double DT= 2.00;
    int Nwr= 20;

    int modYpic[][][]= markpicClrrange(Rycbcr, White, ycbcr, maskrgn);
    ycbcrpic2rgbpic(modYpic, modrgb);



    // +++++++++++++++++++++++++++ Write movie-sequence in random mode ++++++++++++++++++++++++++++

    Writefile=1;
    if(Writefile>0 && Session==0) {

    	String wrfnm= FnameCore+"-"+DfnameCore+"Mx"+Movxy[0]+"My"+Movxy[1]+"Rx"+Rndmv[0];
    	String wname;

    	for(int i=0; i<=Nwr; i++	) {

    		int mskrain[][]= movMaskrgnVec(DT, Movxy, Rndmv, maskrgn);

    		int bufpic[][][]=resample2Dpic(1,1, drsmpic);

    		mapMaskrgn2pic(Ds, mskrain, rsmpic, bufpic);

    		if(i<10) {
    			wname=wrfnm+"R0"+i+".jpg";
    		}
    		else {
    			wname=wrfnm+"R"+i+".jpg";
    		}
    		pictureWriteFile(wname, bufpic);

    		for(int y=0; y< maskrgn.length; y++) {
    			for(int x=0; x< maskrgn[0].length; x++) {
    				maskrgn[y][x]= mskrain[y][x];
    			}
    		}
    	}
    }
    int Mskrain[][]= movMaskrgnVec(DT, Movxy, Rndmv, maskrgn);

    mapMaskrgn2pic(1.0, Mskrain, rsmpic, drsmpic);

     // ++++++++++++++++ Draw pictures ++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Session++;
     super.paintComponent(g);

     img2 = makeImage(rsmpic, W, H);	// paint image of rgbmod to img2
     g2.drawImage(img2, 10, 10, this);	// draw img2 form (10,10)

     img2 = makeImage(modrgb);
     g2.drawImage(img2, 20+rsmpic[0].length, 10, this);

     img2 = makeImage(drsmpic);
     g2.drawImage(img2, 10, 40+rsmpic.length, this);

     g2.drawString("File:"+Fname+" Resmpl-X "+Nx+" Y "+Ny+" Yrange "+Yrange[0]+"--"+Yrange[1]+
    		 		" Yave "+(int)AveMinmax[0][0]+" Ystd "+(int)AveMinmax[0][1], 20,30+H);


  }
}

