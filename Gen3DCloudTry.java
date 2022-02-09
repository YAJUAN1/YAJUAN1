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

public class Gen3DCloudTry extends JPanel
{
  BufferedImage img, imgs;
  Image img2, img3;
  static int Step=0;

  public static void main(String args[])
  {
     JFrame t = new JFrame("Gen3DCloudTry");
     t.getContentPane().add(new Gen3DCloudTry());
     t.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
     t.setSize(2000,1500);
     t.setVisible(true);
  }
  static void plotPrspctv3Dmdl(int clr[], int w, int h, int d, int ld, double thet, int mdl3D[][][], int ppic[][][]) {

	  for(int y=0; y< mdl3D.length; y++) {
		  for(int x=0; x< mdl3D[0].length; x++) {
			  for(int z=0; z< mdl3D[0][0].length; z++) {
				  if(mdl3D[y][x][z]>0) {
					  plotPrspctv(clr, w, h, d, ld, thet, x, y, z, ppic);
				  }
			  }
		  }
	  }
  }
  static void plotPrspctv(int clr[], int w, int h, int d, int ld, double thet, int x, int y, int z, int ppic[][][]) {
	  double lz= ld*((double)z/d);
	  int ix= (int)(x-lz*Math.cos(thet));
	  int iy= (int)(y+lz*Math.sin(thet));
	  if(ix<0 || ix>= ppic[0].length || iy<0 || iy>=ppic.length) {
		  System.out.println("?Illegal perspective (x,y) "+ix+","+iy);
		  return;
	  }
	  movecolor(clr, ppic[iy][ix]);
  }
  static void plotPrspctvSkltn(int clr[], int w, int h, int ld, int dw, int dh, int ppic[][][]) {
	  plotLine(clr, dw,0,     dw,h,      ppic);
	  plotLine(clr, dw,h,     dw+w-1, h, ppic);
	  plotLine(clr, dw+w-1,h, dw+w-1,0,  ppic);
	  plotLine(clr, dw+w-1,0, dw,0,      ppic);

	  plotLine(clr, 0,dh,     0,dh+h-1,  ppic);
	  plotLine(clr, 0,dh+h-1, w,dh+h-1,  ppic);
	  plotLine(clr, w,dh+h-1, w,dh,      ppic);
	  plotLine(clr, w,dh,     0,dh,      ppic);

	  plotLine(clr, dw,0,     0,dh,      ppic);
	  plotLine(clr, dw,h,     0,dh+h-1,  ppic);
	  plotLine(clr, dw+w-1,h, w,dh+h-1,  ppic);
	  plotLine(clr, dw+w-1,0, w,dh,      ppic);
  }
  /* -----------------------------------------------------------------
   * generate XZ picture from 3D model
   *   scl						(x,z)point draw with scl*N(x,z)/H.
   *   clr[c]					showing existence color.
   *   mdl3D[y][x][z]			3D model; 1=Exist, 0= otherwise.
   * OUT
   *   pxz[z][x][c]				XZ picture.
   */
  static void picXZgen( double scl, int clr[], int mdl3D[][][], int pxz[][][]) {

	  for(int x=0; x< mdl3D[0].length; x++) {
		  for(int z=0; z< mdl3D[0][0].length; z++) {
			  double s=0;
			  for(int y=0; y< mdl3D.length; y++) {
				  s+= mdl3D[y][x][z];
			  }
			  movecolor(scl*s/mdl3D.length, clr, pxz[z][x]);
		  }
	  }
  }
  /* --------------------------------------------------------------
   * Setup Z axis existence referring to Y value.
   *   miny, maxy				Prob(p/Depth)= (Y-miny)/(maxy-miny).
   *   rgbpic[y][x][c]			pictuire (rgb based)
   *   vrnd						Y value negative factor. Yval'= Yval - vrnd*random().
   *
   * OUT:
   *   model3D[y][x][z]			3D model for the pic.
   */
  static void pic2Dto3D( double maxy, double miny, double vrnd, int rgbpic[][][], int model3D[][][]) {

	  int val=0;
	  for(int y=0; y<rgbpic.length; y++) {
		  for(int x=0; x< rgbpic[0].length; x++) {
			  val= rgb2y(rgbpic[y][x]);
			  val-= vrnd*Math.random();
			  if(val >= miny)
				  yval2zpoints(miny, maxy, val, model3D[y][x]);
		  }
	  }
  }

  static void pic2Dto3D( double maxy, double miny, int rgbpic[][][], int model3D[][][]) {

	  int val=0;
	  for(int y=0; y<rgbpic.length; y++) {
		  for(int x=0; x< rgbpic[0].length; x++) {
			  val= rgb2y(rgbpic[y][x]);
			  if(val >= miny)
				  yval2zpoints(miny, maxy, val, model3D[y][x]);
		  }
	  }
  }

  static int rgb2y( int c[]) {

		   return (int)(0.29891 * c[0] + 0.58661 * c[1] + 0.11448 * c[2]);
  }
  static void yval2zpoints( double miny, double maxy, int yval, int z[] ) {
	  int n=0;
	  double prb= (yval-miny)/(maxy-miny);
	  for(int k=0; k< z.length; k++) {
		  if(Math.random() <= prb) {
			  z[k]=1;
			  n++;
		  }
	  }
  }
  /* --------------------------------------------------------------
   * Plot an obal on the given pictuire: plotObal()
   *   col[]..				color for plotting
   *   ox, oy..				(ox,oy): point 'O'.
   *   d..					resolution for plotting obal.
   *   a, b..				x**2/a**2 + y**2/b**2 = 1.
   *
   * OUT:
   *   pic[][][]
   */
  static void plotObal(int col[], int ox, int oy, double d, double a, double b, int pic[][][]) {

	  double dy[]= {0,0};
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);

		  int ix,iy;
		  ix=(int)(ox+dx);
		  iy=(int)(oy-dy[0]);

		  if(iy<0 || iy>=pic.length || ix<0 || ix>= pic[0].length) 	continue;
		  movecolor(col, pic[iy][ix]);

		  iy=(int)(oy-dy[1]);
		  if(iy<0 || iy>=pic.length) 		continue;
		  movecolor(col, pic[iy][ix]);
	  }
  }
  /* -------------------------------------------------------------------------------
   * Plot normal vectors of obal: plotObalNormalVecLen().
   *   col[]..						color.
   *   ox,oy..						point O for the obal
   *   d..							plot interval (x-axis).
   *   inf..						dydx max. value.
   *   len..						vector len.
   *   a,b..						x^2/a^2+y^2/b^2=1.
   * OUT:
   *   pic[y][x][c]
   *
   */
  static void plotObalNormalVecLen(int col[], int ox, int oy, double d, double inf, double len,
		  							double a, double b, int pic[][][]) {

	  double dy[]= {0,0}, dydx=0 ;
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);
		  dydx= fobalNormalvec( dx, a, b);

		  plotVecLen(col, (int)(ox+dx), (int)(oy-dy[0]),  dydx, inf, len, pic);

		  plotVecLen(col, (int)(ox+dx), (int)(oy+dy[0]), -dydx, inf, len, pic);
	  }
  }
  static void plotObalNormalVecRndm(int col[], int ox, int oy, double d, double inf, double len,
									double a, double b, int pic[][][]) {

	  double dy[]= {0,0}, dydx=0 ;
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);
		  dydx= fobalNormalvec( dx, a, b);

		  plotVecRndmLen(col, (int)(ox+dx), (int)(oy-dy[0]),  dydx, inf, len, pic);

		  plotVecRndmLen(col, (int)(ox+dx), (int)(oy+dy[0]), -dydx, inf, len, pic);
	  }
  }
  /* -------------------------------------------------------------------------------
   * Plot normal vectors of obal in triangle-form: plotObalNormalVecTangl().
   *   col[]..						color.
   *   w..							width of triangle.
   *   r..							shrink rate for traiangle shape.
   *   ox,oy..						point O for the obal.
   *   d..							plot interval (x-axis).
   *   inf..						dydx max. value.
   *   len..						vector len.
   *   a,b..						x^2/a^2+y^2/b^2=1.
   * OUT:
   *   pic[y][x][c]
   */
  static void plotObalNormalVecTangl(int col[], double cvar[], int coln[], int w, double r, int ox, int oy, double d, double inf,
										double len, double a, double b, int pic[][][]) {

	  double dy[]= {0,0}, dydx;
	  int c[]= new int[col.length];

	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);						// get Y values in dy[0/1].

		  dydx= fobalNormalvec( dx, a, b);			// compute normal vector(slope) at dx.
		  int s= (dx<0)? 1:0 ;

		  double rc= Math.max( cvar[0], cvar[1]*Math.random());
		  movecolor( rc, col, c);

		  plotVecRndmTangl(c, coln, s, w, r, (int)(ox+dx), (int)(oy-dy[0]),  dydx, inf, len, pic);

		  plotVecRndmTangl(c, coln, s, w, r, (int)(ox+dx), (int)(oy+dy[0]), -dydx, inf, len, pic);
	  }
	}

  static void plotObalNormalVecTangl(int col[], int w, double r, int ox, int oy, double d, double inf, double len,
									double a, double b, int pic[][][]) {

	  double dy[]= {0,0}, dydx;
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);						// get Y values in dy[0/1].

		  dydx= fobalNormalvec( dx, a, b);			// compute normal vector(slope) at dx.

		  plotVecRndmTangl(col, w, r, (int)(ox+dx), (int)(oy-dy[0]),  dydx, inf, len, pic);

		  plotVecRndmTangl(col, w, r, (int)(ox+dx), (int)(oy+dy[0]), -dydx, inf, len, pic);
	  }
  }
  static void plotObalNormalVecTangl(int col[], int coln[], int w, double r, int ox, int oy, double d, double inf,
		  							double len, double a, double b, int pic[][][]) {

	  double dy[]= {0,0}, dydx;
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);						// get Y values in dy[0/1].

		  dydx= fobalNormalvec( dx, a, b);			// compute normal vector(slope) at dx.
		  int s= (dx<0)? 1:0 ;

		  plotVecRndmTangl(col, coln, s, w, r, (int)(ox+dx), (int)(oy-dy[0]),  dydx, inf, len, pic);

		  plotVecRndmTangl(col, coln, s, w, r, (int)(ox+dx), (int)(oy+dy[0]), -dydx, inf, len, pic);
	  	}
  }

  static void plotObalNormalVec(int col[], int ox, int oy, double d, double lx, double a, double b, int pic[][][]) {

	  double dy[]= {0,0}, dydx=0, ly=0, llx=0 ;
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);
		  dydx= fobalNormalvec( dx, a, b);
		  llx=lx;
		  if(dx<0) {
			  llx=-lx;
		  }
		  ly= dydx*llx;
		  //					System.out.println(dx+","+dy[0]+" df "+dydx+" llx "+llx+" ly "+ly);
		  int ix,iy, hx, hy;
		  ix=(int)(ox+dx);
		  iy=(int)(oy-dy[0]);
		  if(dx == 0) {
			  hx=ix;
			  hy=(int)(iy-lx);
		  }
		  else {
			  hx=(int)(ox+dx+llx);
			  hy=(int)(oy-dy[0]-ly);
		  }
		  if(iy<0 || iy>=pic.length || ix<0 || ix>= pic[0].length) 	continue;
		  if(hy<0 || hy>=pic.length || hx<0 || hx>= pic[0].length)	continue;

		  plotLine(col, ix, iy, hx, hy, pic);

		  iy=(int)(oy-dy[1]);
		  dydx=-dydx;
		  ly= -ly;
		  hy=(int)(oy-dy[1]-ly);

		  if(iy<0 || iy>=pic.length || hy<0 || hy>=pic.length) 		continue;

		  plotLine(col, ix, iy, hx, hy, pic);
	  }
  }
  /*---------------------------------------------------------------------------
   * Plot a vector of random length: plotVecRndmLen()
   *   col[]..					color vector.
   *   ox,oy..					starting point (x,y).
   *   dydx..					tangent of the vector.
   *   inf..					infinity number ex. 9999.
   *   len..					vector length.
   * OUT:
   *   pic[][][]..				picture campus to be drawn.
   */
  static void plotVecRndmLen(int col[], int ox, int oy, double dydx, double inf, double len, int pic[][][]) {

	  double rl= Math.random()*len - len/2.0 ;
	  int s = (rl<0)? -1:1;

	  int ix1, ix2, iy1, iy2;
	  double dx;
	  if(Math.abs(dydx) >= inf) {
		  ix1=ox;
		  ix2=ox;
		  iy1=oy;
		  iy2=(int)(oy+rl);
	  }
	  else if (dydx == 0){
		  ix1=ox ;
		  ix2=(int)(ox+rl);
		  iy1=oy;
		  iy2=oy;
	  }
	  else {
		  dx= rl/Math.sqrt(1+dydx*dydx);
		  ix1= ox ;
		  ix2=(int)(ox+dx);
		  iy1= oy ;
		  iy2=(int)(oy-dx*dydx);
	  }
	  plotLine(col, ix1, iy1, ix2, iy2, pic);
  }
  /* -------------------------------------------------------------------------------
   * Plot normal vector in triangle form: plotVecRndmTangl().
   *   col[]..						color.
   *   w..							width of triangle.
   *   r..							shrink rate for traiangle shape.
   *   rl..							length of the vector.
   *   ox,oy..						point O for the obal.
   *   d..							plot interval (x-axis).
   *   inf..						dydx max. value.
   *   len..						vector len.
   * OUT:
   *   pic[y][x][c]
   */
  static void plotVecRndmTangl(int col[], int w, double r, int ox, int oy, double dydx, double inf, double len,
		  						int pic[][][]) {

	  double rl= Math.random()*len - len/2.0 ;

	  if(Math.abs(dydx) >= inf) {						// ..... Vertical case.
		  plotVecTanglVH(0, col, ox, oy, w, r, rl, pic);
	  }
	  else if (dydx == 0){								// ..... Horizontal case.
		  plotVecTanglVH(1, col, ox, oy, w, r, rl, pic);
	  }
	  else {
		  plotVecTangl(col, ox, oy, dydx, w, r, rl, pic);
	  }
  }
  static void plotVecRndmTangl(int col[], int coln[], int s, int w, double r, int ox, int oy, double dydx, double inf,
		  						double len, int pic[][][]) {

	  double rl= Math.random()*len - len/2.0 ;

	  int cl[]= {0,0,0};

	  if(s==1 ) {
		  if(rl < 0)			movecolor(col, cl);
		  else					movecolor(coln, cl);
	  }
	  else {
		  if(rl < 0)			movecolor(coln, cl);
		  else					movecolor(col, cl);
	  }

	  if(Math.abs(dydx) >= inf) {						// ..... Vertical case.
		  plotVecTanglVH(0, cl, ox, oy, w, r, rl, pic);
	  }
	  else if (dydx == 0){								// ..... Horizontal case.
		  plotVecTanglVH(1, cl, ox, oy, w, r, rl, pic);
	  }
	  else {
		  plotVecTangl(cl, ox, oy, dydx, w, r, rl, pic);
	  }
  }
  static void plotVecTangl( int col[], int ox, int oy, double dydx, int w, double r, double rl, int pic[][][]) {

	  double dx= rl/Math.sqrt(1+dydx*dydx);
	  													// normal vector plotting (the center of the triangle.)
	  plotLine(col, ox, oy, (int)(ox+dx), (int)(oy-dx*dydx), pic);

	  double p= -1.0/dydx;			// p.. tangent vector

	  if(Math.abs(p) <= 1.0) {		// compute triangle normal vectors in X-first manner.
		  for(double i=0.5 ; i<= w; i+=0.5) {

			  rl*=r;
			  dx= rl/Math.sqrt(1+dydx*dydx);

			  plotLine( col, (int)(ox+i), (int)(oy-p*i), (int)(ox+i+dx), (int)(oy-p*i-dx*dydx), pic);
			  plotLine( col, (int)(ox-i), (int)(oy+p*i), (int)(ox-i+dx), (int)(oy+p*i-dx*dydx), pic);
		  }
	  }
	  else {						// compute triangle normal vectors in Y-first manner.
		  p= 1.0/p;
		  double dxdy=1./dydx,  dy;
		  for(double k=0.5; k<= w; k+=0.5) {

			  rl*=r;
			  dy= rl/Math.sqrt(1+dxdy*dxdy);
			  if(dxdy < 0) dy= -dy;

			  plotLine( col, (int)(ox+p*k), (int)(oy-k), (int)(ox+p*k+dy*dxdy), (int)(oy-k-dy), pic);
			  plotLine( col, (int)(ox-p*k), (int)(oy+k), (int)(ox-p*k+dy*dxdy), (int)(oy+k-dy), pic);

		  }

	  }

  }
  /* -------------------------------------------------------------------------------
   * Plot normal vector in triangle form w/ pos/nega colors: plotVecRndmTangl().
   *   col[],coln[]..				colors for positive/negative vectors.
   *   w..							width of triangle.
   *   r..							shrink rate for traiangle shape.
   *   rl..							length of the vector.
   *   ox,oy..						point O for the obal.
   *   d..							plot interval (x-axis).
   *   inf..						dydx max. value.
   *   len..						vector len.
   * OUT:
   *   pic[y][x][c]
   */
  static void plotVecTangl( int col[], int coln[], int ox, int oy, double dydx, int w, double r, double rl, int pic[][][]) {

	  double dx= rl/Math.sqrt(1+dydx*dydx);
	  													// normal vector plotting (the center of the triangle.)
	  plotLine(col, ox, oy, (int)(ox+dx), (int)(oy-dx*dydx), pic);

	  double p= -1.0/dydx;			// p.. tangent vector

	  if(Math.abs(p) <= 2.0) {		// compute triangle normal vectors in X-first manner.
		  for(double i=0.5 ; i<= w; i+=0.5) {

			  rl*=r;
			  dx= rl/Math.sqrt(1+dydx*dydx);

			  plotLine( col, (int)(ox+i), (int)(oy-p*i), (int)(ox+i+dx), (int)(oy-p*i-dx*dydx), pic);
			  plotLine( col, (int)(ox-i), (int)(oy+p*i), (int)(ox-i+dx), (int)(oy+p*i-dx*dydx), pic);
		  }
	  }
	  else {						// compute triangle normal vectors in Y-first manner.
		  p= 1.0/p;
		  double dxdy=1./dydx,  dy;
		  for(double k=0.5; k<= w; k+=0.5) {

			  rl*=r;
			  dy= rl/Math.sqrt(1+dxdy*dxdy);
			  if(dxdy < 0) dy= -dy;

			  plotLine( col, (int)(ox+p*k), (int)(oy-k), (int)(ox+p*k+dy*dxdy), (int)(oy-k-dy), pic);
			  plotLine( col, (int)(ox-p*k), (int)(oy+k), (int)(ox-p*k+dy*dxdy), (int)(oy+k-dy), pic);

		  }

	  }

  }
  static void plotVecTanglVH( int vh, int col[], int ox, int oy, int w, double r, double rl, int pic[][][]) {

	  int ix1, ix2, iy1, iy2;
	  if(vh==0) {	// Vertical Lines

		  plotLine(col, ox, oy, ox, (int)(oy-rl), pic);
		  for(int i=1; i<=w; i++) {
			  rl*=r;
			  plotLine(col, ox+i, oy, ox+i, (int)(oy-rl), pic);
			  plotLine(col, ox-i, oy, ox-i, (int)(oy-rl), pic);
		  }
	  }
	  else {		// Horizontal Lines

		  plotLine(col, ox, oy, (int)(ox+rl), oy, pic);
		  for(int i=1; i<=w; i++) {
			  rl*=r;
			  plotLine(col, ox, oy-i, (int)(ox+rl), oy, pic);
			  plotLine(col, ox, oy+i, (int)(ox+rl), oy, pic);
		  }
	  }
  }
  /*---------------------------------------------------------------------------
   * Plot a vector of given length: plotVecLen()
   *   col[]..					color vector.
   *   ox,oy..					starting point (x,y).
   *   dydx..					tangent of the vector.
   *   inf..					infinity number ex. 9999.
   *   len..					vector length.
   * OUT:
   *   pic[][][]..				picture campus to be drawn.
   */
  static void plotVecLen(int col[], int ox, int oy, double dydx, double inf, double len, int pic[][][]) {

	  int ix1, ix2, iy1, iy2;
	  double dx;
	  if(Math.abs(dydx) >= inf) {
		  ix1=ox;
		  ix2=ox;
		  iy1=(int)(oy-len);
		  iy2=(int)(oy+len);
	  }
	  else if (dydx == 0){
		  ix1=(int)(ox-len);
		  ix2=(int)(ox+len);
		  iy1=oy;
		  iy2=oy;
	  }
	  else {
		  dx= len/Math.sqrt(1+dydx*dydx);
		  ix1=(int)(ox-dx);
		  ix2=(int)(ox+dx);
		  iy1=(int)(oy+dx*dydx);
		  iy2=(int)(oy-dx*dydx);
	  }
	  plotLine(col, ix1, iy1, ix2, iy2, pic);
  }
  /*
   * draw & paint obal with a given color.
   *   col[c]				color max values.
   *   cvar[rmin/rmax]		final color = col[c]x(Max(cvar[0], cvar[1] x random())).
   */
  static void plotObalFill(int col[], double cvar[], int ox, int oy, double d, double a, double b, int pic[][][]) {

	  double dy[]= {0,0};
	  int c[]= new int [col.length];
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);

		  int ix,iy;
		  ix=(int)(ox+dx);
		  iy=(int)(dy[0]+oy);

		  if(iy<0 || iy>=pic.length || ix<0 || ix>= pic[0].length) 	continue;

		  int iy2 =(int)(dy[1]+oy);
		  if(iy2<0 || iy2>=pic.length) 		continue;
		  double var= Math.max(cvar[0], cvar[1]*Math.random());
		  movecolor(var, col, c);
		  plotLine( c, ix, iy, ix, iy2, pic);
	  }
  }
  static void plotObalFill(int col[], int ox, int oy, double d, double a, double b, int pic[][][]) {

	  double dy[]= {0,0};
	  for( double dx = -a; dx <= a; dx+=d) {

		  fobal( dx, a, b, dy);


		  int ix,iy;
		  ix=(int)(ox+dx);
		  iy=(int)(dy[0]+oy);

		  if(iy<0 || iy>=pic.length || ix<0 || ix>= pic[0].length) 	continue;

		  int iy2 =(int)(dy[1]+oy);
		  if(iy2<0 || iy2>=pic.length) 		continue;

		  plotLine( col, ix, iy, ix, iy2, pic);
	  }
  }
  /* --------------------------------------------------------------
   * function obal: fobal()
   *   x..			x.
   *   a,b..		x**2/a**2 + y**2/b**2 = 1.
   * OUT:
   *   dy[0/1]..	y= +- b sqrt[ x**2/a**2 - 1 ].
   */
  static void fobal(double x, double a, double b, double y[]	) {

	  y[0]= b*Math.sqrt(1-x*x/(a*a));
	  y[1]= -y[0];

  }
  static double fobalNormalvec(double x, double a, double b) {

	  return a*a*Math.sqrt(1-x*x/(a*a))/(b*x);


  }
  static void plotLine(int col[], int x1, int y1, int x2, int y2, int rgb[][][]) {

	  double dy, dx, rmy, rmx;
	  if (x1<0 || y1<0 || x2<0 || y2<0) {
		  System.out.println(" Index < 0 Error.");
		  return;
	  }
	  if (x1>=rgb[0].length || x2>=rgb[0].length || y1>=rgb.length || y2>=rgb.length) {
		  System.out.println(" Index exceeds max.");
		  return;
	  }

	  dy = y2-y1;   dx = x2-x1;
	  if (dy==0){
		  if(x1 > x2){
			  for(double j=x2; j<=x1; j+=0.5){

				  plotDot(col, (int)(j+0.5), y1, rgb);
			  }
		  }
		  else {
			  for (int j=x1; j<=x2; j++){
				  plotDot(col, j, y1, rgb);
			  }
		  }
		  return;
	  }
	  else if (dx==0){
		  if(y1 > y2){
			  for(double i=y2; i<=y1; i+=0.5){
				  plotDot(col, x1, (int)(i+0.5), rgb);
			  }
		  }
		  else {
			  for(double i=y1; i<=y2; i+=0.5){
				  plotDot(col, x1, (int)(i+0.5), rgb);
			  }
		  }
		  return;
	  }
	  rmy = dy/dx;
	  rmx = dx/dy;
	  int lx,ly;
	  if(Math.abs(dy)<Math.abs(dx)){
		  if(dx < 0){
			  for(double i=dx; i<=0; i+=0.5){
				  lx=(int)(x1+i);
				  ly=(int)(rmy*i+y1);
				  plotDot(col, lx, ly, rgb);
			  }
		  }
		  else {
			  for(double i=0; i<=dx; i+=0.5){
				  lx=(int)(x1+i);
				  ly=(int)(rmy*i+y1);
				  plotDot(col, lx, ly, rgb);
			  }
		  }
		  return;
	  }
	  if (dy < 0){
		  for (double i=dy; i<=0; i+=0.5){
			  lx=(int)(rmx*i+x1);
			  ly=(int)(y1+i);
			  plotDot(col, lx, ly, rgb);
		  }
	  }
	  else {
		  for (double i=0; i<= dy; i+=0.5){

			  lx=(int)(rmx*i+x1);
			  ly=(int)(y1+i);
			  plotDot(col, lx, ly, rgb);
		  }
	  }
  }
  public static void plotPlus(int col[], int x, int y, int rgb[][][]){
	  int d=3;
	  plotLine(col, x-d, y, x+d, y, rgb);
	  plotLine(col, x, y-d, x, y+d, rgb);
  }
  public static void overpaintColor(int clr[], double wc, int px[]){
	  for(int c=0; c<px.length; c++){
		  px[c]= (int)(wc*clr[c]+(1-wc)*px[c]+0.5);
		  px[c]= Math.min(255, Math.max(0, px[c]));
	  }
  }
  public static void plotDot(int col[], int x, int y, int rgb[][][]){
	  if(x<0 || y<0 || x>=rgb[0].length || y>=rgb.length)   return;
	  copyRgb(col, rgb[y][x]);
  }
  public static void copyRgb(int p1[], int p2[]){
	  for (int c=0; c< 3; c++) p2[c]= p1[c];
  }

  static void applyOperator2Pic(double op[][], int pics[][][], int picd[][][]){
	  int dy= op.length/2;
	  int dx= op[0].length/2;

	  for(int y=dy; y< pics.length-dy; y++){
		  for(int x=dx; x< pics[0].length-dx; x++){
			  for(int clr=0; clr< pics[0][0].length; clr++){
				  double v = applyOperator(op, x, y, clr, pics);
				  picd[y][x][clr]= (int)v;
			  }
		  }
	  }
  }
  static double applyOperator(double op[][], int x, int y, int c, int pic[][][]){
	  double v=0;
	  int cy = op.length/2;
	  int cx = op[0].length/2;

	  for(int i=0; i< op.length; i++){
		  int dy= y+ (i-cy);
		  for(int j=0; j< op[0].length; j++){
			  int dx= x+ (j-cx);
			  v += op[i][j]*pic[dy][dx][c];
		  }
	  }
	  return v;
  }
  	static void fillcolor2Pic(int clr[], int pic[][][]){
  		for(int y=0; y< pic.length; y++){
  			for(int x=0; x< pic[0].length; x++){
  				movecolor(clr, pic[y][x]);
  			}
  		}
  	}
	/* -----------------------------------------------------
	 * Gaussian Filter Generation: gaussmat()
	 *   l				l x l matrix filter
	 *   var			variance
	 * RETRUN
	 *   mat[][]		gaussian fileter of l x l matrix
	 */
    static double[][] gaussmat(int l, double var){
    	double mat[][]= new double[l][l];
		int d=l/2;
		if (mat.length!=(d*2+1) || mat[0].length!=(d*2+1)){
			System.out.println(" Gaussian matrix error:"+mat.length);
			return null;
		}
		for(int x=-d; x<=d; x++){
			for(int y=-d; y<=d; y++){
				mat[y+d][x+d]= Math.exp(-(x*x+y*y)/(2*var))/(2*Math.PI*var);
			}
		}
		double a=0;
		for(int i=0; i<mat.length; i++){
			for(int j=0; j<mat.length; j++){
				System.out.print(" "+(int)(100*mat[i][j]));
				a+=mat[i][j];
			}
			System.out.println();
		}
		System.out.println("*sum:"+a);
		for(int y=0; y< mat.length; y++){
			for(int x=0; x< mat[0].length; x++){
				mat[y][x]/= a;
			}
		}
		return mat;
    }
	static void gaussmat(int l, double var, double mat[][]){
		int d=l/2;
		if (mat.length!=(d*2+1) || mat[0].length!=(d*2+1)){
			System.out.println(" Gaussian matrix error:"+mat.length);
			return;
		}
		for(int x=-d; x<=d; x++){
			for(int y=-d; y<=d; y++){
				mat[y+d][x+d]= Math.exp(-(x*x+y*y)/(2*var))/(2*Math.PI*var);
			}
		}
		double a=0;
		for(int i=0; i<mat.length; i++){
			for(int j=0; j<mat.length; j++){
				System.out.print(" "+(int)(100*mat[i][j]));
				a+=mat[i][j];
			}
			System.out.println();
		}
		System.out.println("*sum:"+a);
	}

	int [][][] applyOperator2Pic(double oprt[][], int pic[][][]){

		int pout[][][] = new int[pic.length][pic[0].length][pic[0][0].length] ;
		int loprt= oprt.length;
		int dl = loprt/2;
		for(int y=dl; y< pic.length-dl; y++){
			for(int x=dl; x< pic[0].length-dl; x++){
				applyOperator1(x, y, oprt, pic, pout);
			}
		}
		return pout;
	}
	int [][] applyOperator2Pic(double oprt[][], int pic[][]){
		int pout[][] = new int[pic.length][pic[0].length] ;
		int loprt= oprt.length;
		int dl = loprt/2;
		for(int y=dl; y< pic.length-dl; y++){
			for(int x=dl; x< pic[0].length-dl; x++){
				applyOperator1(x, y, oprt, pic, pout);
			}
		}
		return pout;
	}
	void applyOperator1(int x0, int y0, double oprt[][], int pic[][][], int pout[][][]){
		int dl= oprt.length/2;
		for(int c=0; c< pout[0][0].length; c++) {
			double v=0;
			for(int y= y0-dl; y <= y0+dl; y++){
				int oy= y-(y0-dl);
				for(int x= x0-dl; x <= x0+dl; x++){
					int ox= x-(x0-dl);
					v+=oprt[oy][ox]*pic[y][x][c] ;
				}
			}
			pout[y0][x0][c] = (int)v;
		}
	}
	void applyOperator1(int x0, int y0, double oprt[][], int pic[][], int pout[][]){
		int dl= oprt.length/2;
		int v=0;
		for(int y= y0-dl; y <= y0+dl; y++){
			int oy= y-(y0-dl);
			for(int x= x0-dl; x <= x0+dl; x++){
				int ox= x-(x0-dl);
				v+=(int)(oprt[oy][ox]*pic[y][x]);
			}
		}
		pout[y0][x0] = v;
	}
	static void showmat(int mat[][]){
		for(int y=0; y<mat.length; y++){
			for(int x=0; x<mat[0].length; x++){
				System.out.printf(" %3d", mat[y][x]);
			}
			System.out.println();
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
    			movevec(pic[iy][ix], ps[y][x]);
    			x++;
    		}
    		y++;
    	}
    	return ps;
    }
    static void movevec(int v1[], int v2[]){
    	for(int i=0; i< v1.length; i++) 	v2[i]= v1[i];
    }
    static void ycbcrpic2rgbpic(int ycbcr[][][], int rgb[][][]){
    	for(int y=0; y<ycbcr.length; y++){
    		for(int x=0; x<ycbcr[0].length; x++){
    			ycbcr2rgb(ycbcr[y][x], rgb[y][x]) ;
    		}
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
    static int [][][] grayscalepic2rgbpic(int gsc[][]){
    	int rgbp[][][]= new int [gsc.length][gsc[0].length][3];
    	for(int y=0; y<gsc.length; y++){
    		for(int x=0; x<gsc[0].length; x++){
    			for(int c=0; c< rgbp[0][0].length; c++){
    				rgbp[y][x][c]= gsc[y][x];
    			}
    		}
    	}
    	return rgbp;
    }
    static int[][] rgbpic2grayscale(int rgb[][][]){
    	int gs[][]= new int[rgb.length][rgb[0].length];

    	for(int y=0; y<rgb.length; y++){
    		for(int x=0; x<rgb[0].length; x++){
    			gs[y][x]= rgb2grayscale(rgb[y][x]) ;
    		}
    	}
    	return gs;
    }
    static int rgb2grayscale(int rgb[]){
    	return (int)(0.29891 * rgb[0] + 0.58661 * rgb[1] + 0.11448 * rgb[2]);
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
	Image makeGrayscaleImage(int gspic[][]) {
		int h= gspic.length;
		int w= gspic[0].length;
		int [] pixel = new int[w*h];

		for(int i=0; i<h;i++){
			for(int j=0; j<w; j++) {
				pixel[j+i*w]= (255<<24) | (gspic[i][j]<<16) |
					(gspic[i][j]<<8) | (gspic[i][j]);
			}
		}
		return(createImage(new MemoryImageSource(w, h, pixel, 0,w)));
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
	static void movecolor( double r, int clr[], int p[]) {
		for(int c=0; c< p.length; c++) {
			p[c]=(int)(r*clr[c]);
			p[c]= Math.max(0, Math.min(255, p[c]));
		}
	}
	static void movevec(double c1[], int c2[]){
		for(int c=0; c<c1.length; c++){
			c2[c]= (int)c1[c];
		}
	}

  public void update(Graphics g)
  {
    paint(g);
  }

  /***********************************************************
   * Cloud Generation Project: plot obal
   *
   ***********************************************************/
  public void paintComponent(Graphics g)
  {
    Graphics2D g2= (Graphics2D)g;

    int Black[]	= {  0,   0,   0};			// Edge drawing color.
    int White[] = {255, 255, 255};			// background color.
    int Red[]  	= {255,   0,   0};
    int Green[]	= {   0,255,   0};

    String FnameCore="海岸03";
    String Fname=FnameCore+".jpg";		// Input picture file name
    int Nx= 1, Ny=Nx;					// Resampling rate.

    boolean WriteFile= true;				// write file or not

    //------------------------ Picture Parameters -------------------------------
    int W = 		 400;					// Width
    int H =	    	 300;					// Height
    int D =			 150;					// Depth
    int LD=			 250;					// Lz in sloped Depth
    int XO =		W/2;					// (x,y) for "O".
    int YO =		H/2;
    double Ymax=    230;					// Y max value.
    double Ymin=     50;					// Y min. value as clouds.
    double Zscl=	 15;
    double Thet=	 Math.toRadians(35);	// sloped DL angle to horizontal X axis.
    int DY= (int)(Math.sin(Thet)*LD);		// Extra Hight for sloped Z in perspective picture.
    int DX= (int)(Math.cos(Thet)*LD);		// Extra Width for sloped Z in perspective picture.

    double Vrnd=        255/3;				// Suppress factor for Y value when creating 3D model.
    double Vcol[]= {0.2, 1 };				// Rmin, Rmax for painting color.

    //------------------------ Obal Parameters & Plot Obal ---------------------
    double Aobl =    160;					// x^2/A^2 + y^2/B^2 = 1.
    double Bobl =	 100;
    double Rsl  =	 0.2;					// resolution for drawing obal.
    double Dx   =      5;					// vector plot interval in X.
    double Len  =     60;					// vector length.
    double Inf  =    999;					// infinity number.
    int    WL   =      13;					// width of triangle lines.
    double RL   =   0.85;					// triangle shrink rate.

    int cldpic[][][]= 	new int[H][W][3];			// picture for cloud drawing
    int cld3D[][][]=	new int[H][W][D];			// cloud 3D model
    int prspic[][][]=	new int[H+DY+1][W+DX+1][3];		// perspective cloud picture

    plotObalFill(White, Vcol, XO, YO, Rsl, Aobl, Bobl, cldpic);										// draw Obal.
    plotObalNormalVecTangl(White, Vcol, Black, WL, RL, XO, YO, Dx, Inf, Len, Aobl, Bobl, cldpic);		// draw Random vec. triangle.

    int Lgauss     =  7 ;				// gaussian window length.
    double Vgauss  =  8.00;				// gaussian window variance.
    int Ngauss     = 55;				// gaussian apply loop.

    double Gauflt[][] = gaussmat(Lgauss, Vgauss);

  //------------------------ Gaussian Fileter & Generate 3D Cloud ---------------------
    int gspic[][][]= applyOperator2Pic(Gauflt, cldpic);
    for(int k=0; k<Ngauss; k++)
    	gspic= applyOperator2Pic(Gauflt, gspic);

    pic2Dto3D(Ymax, Ymin, Vrnd, gspic, cld3D);
    int uppic[][][]= new int[D][W][3];			// upper XZ picture.
    picXZgen(Zscl, White, cld3D, uppic);

    plotPrspctvSkltn(Green, W, H, LD, DX, DY, prspic);
    plotPrspctv3Dmdl(White, W, H, D, LD, Thet, cld3D, prspic);

    if(WriteFile && Step==0) {
    	String wrname="Obl"+W+"x"+H+"x"+D+"Th"+(int)Math.toDegrees(Thet)+"A"+(int)Aobl+"B"+(int)Bobl+"L"+(int)Len+"WL"+WL+
    			"NG"+Ngauss+"Vmin"+Vcol[0]+"Vmax"+Vcol[1]+"Vsp"+(int)Vrnd;
    	pictureWriteFile(wrname+"P.jpg", prspic);
    	pictureWriteFile(wrname+"F.jpg", gspic);
    	pictureWriteFile(wrname+"S.jpg", cldpic);
    }

    //int spic[][][]    = readPicrgb(Fname);				// Read a input picture file Fname
    //int tpic[][][]	  = resample2Dpic(Ny, Nx, spic);

    //-------------------------- Drawing Pictures -------------------------------
     super.paintComponent(g);
     Step++;

     img2 = makeImage(cldpic);
     g2.drawImage(img2, 10, 10, this);

     img2 = makeImage(gspic);
     g2.drawImage(img2, 20+W, 10, this);

     img2 = makeImage(uppic);
     g2.drawImage(img2, 10, 20+H, this);

     img2 = makeImage(prspic);
     g2.drawImage(img2, 20+W, 20+H, this);

     g2.drawString(" WxHxD "+W+"x"+H+"x"+D+" Ld "+LD+" THT "+Math.toDegrees(Thet)+" A,B "+Aobl+","+Bobl, 10,40+H+D);
     g2.drawString(" Vlen "+Len+" wL "+WL+" Rmin,Rmax "+Vcol[0]+","+Vcol[1]+" Vsup "+(int)Vrnd, 10,60+H+D);
  }
}

