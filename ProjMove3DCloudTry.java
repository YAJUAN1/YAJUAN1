package cloudRain;

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

public class ProjMove3DCloudTry extends JPanel
{
  BufferedImage img, imgs;
  Image img2, img3;
  static int Step=0;

  public static void main(String args[])
  {
     JFrame t = new JFrame("ProjMove3DCloudTry");
     t.getContentPane().add(new ProjMove3DCloudTry());
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

	  double dx= ld*Math.cos(thet);
	  double lz= ld*((double)z/d);

	  int ix= (int)(dx+x-lz*Math.cos(thet));
	  int iy= (int)(y+lz*Math.sin(thet));

	  if(ix<0 || ix>= ppic[0].length || iy<0 || iy>=ppic.length) {
		  System.out.println("?Illegal perspective (x,y) "+ix+","+iy);
		  return;
	  }
	  movecolor(clr, ppic[iy][ix]);
  }
  /*
   * Projection 3D cloud model to a 2D picture
   *   clr[]						color vector to draw.
   *   rc							resolve factor.
   *   ld							z axis vertual length
   *   thet							z axis angle to horizontal line.
   *   mdl3D[y][x][z]				3D model to be mapped.
   *   x0/y0, x1/y1					relative box for projection in ppic.
   * OUT
   *   ppic[y][x][c]				projection picture.
   */
  static void prjct3Dmdl2pic(int clr[], double rc, int ld, double thet,
		  int mdl3D[][][], double x0, double y0, double x1, double y1, int ppic[][][]) {

	  for(int y=0; y< mdl3D.length; y++) {
		  for(int x=0; x< mdl3D[0].length; x++) {
			  for(int z=0; z< mdl3D[0][0].length; z++) {
				  if(mdl3D[y][x][z]>0) {
					  double r= Math.random()*rc;
					  prjct3Dpoint2pic(r, clr, mdl3D[0].length, mdl3D.length, mdl3D[0][0].length,
							  ld, thet, x, y, z, x0,y0, x1,y1,ppic);
				  }
			  }
		  }
	  }
  }
  /*
   * Projection 3D cloud model to a 2D picture.
   *   clr[]					cloud color.
   *   w/h/d					width/height/depth for 3D model
   *   ld						depth length in pixel.
   *   thet						z angle to the horizontal line.
   *   x/y/z					the point to be mapped on to the picture.
   *   x0,y0					the projection box left-up most point in relative axis(0-1.0)
   *   x1,y1					the projection box right-down most point.
   * OUT
   *   p[y][x][c]				the 2D picture for output.
   *
   *
   */
  static void prjct3Dpoint2pic(double cr, int clr[], int w, int h, int d, int ld, double thet, int x, int y, int z,
		 double x0, double y0, double x1, double y1, int p[][][]) {

	  double dx = ld*Math.cos(thet);
	  double dy = ld*Math.sin(thet);
	  double aw = dx+w;
	  double ah = dy+h;

	  double lz= ld*((double)z/d);

	  int ix= (int)(dx+x-lz*Math.cos(thet));
	  int iy= (int)(y+lz*Math.sin(thet));

	  double rx= ix/aw;
	  double ry= iy/ah;

	  double pw= p[0].length;
	  double ph= p.length;

	  int ix0= (int)(x0*pw);
	  int iy0= (int)(y0*ph);
	  int ix1= (int)(x1*pw);
	  int iy1= (int)(y1*ph);

	  int px= (int)(ix0 + (ix1-ix0)*rx);
	  int py= (int)(iy0 + (iy1-iy0)*ry);

	  if(px<0 || px>=pw || py<0 || py>= ph) {
		  System.out.println("?Illegal projction x,y: "+px+","+py);
		  return;
	  }
	  movecolor(cr, clr, p[py][px]);
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
	  double prb= (yval-miny)/maxy; 	// prev. code: (yval-miny)/(maxy-miny)

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
  static void plotVecTangl( int col[], int coln[], int o��&�	U4A�N�q"����g�X��57�/�tM*IxE���q��?R	:OL��!��n��-E�O5�W�����Y��noe��ˉ�-���Q=��-I�,8��;���夋���m�����N��z�?�W�wi�YM��* Fˡ�~&)P��N2^7�)׏�����`Zb	�`H3���.~�x�Vˋ{�.\MT���$��.�i3�R^*��I�W>��(��>�W]l�5����E?^��͔g$a�,	��߶X3[Z�5������Ĭ�G=ួPn~x��V6����a�y�H�fR�czP���D�]c��A�"��g��ǟC���<w�+v�s-;���<D�D2�:#n;�%���U80N/p �j^db�ƶ����h_����UL�k]�$g]�M���Y�+t;���>+yvR1�}�5�4�Q�:!)S�7���C��g�
"(�S9r����_
�9��X�[����3%���*�,*�"���䇃�jĲP2uL�2�(b��� p�BBg�k�T���
��,�>F�+4���
�8�rߖ�T"aS��R��'e�1�bs!�z='9"�e?���wU�����v�w�<��{oId(=�z�0ڔ�P�ą�A�����2�ͥ.Eͭ^�}�l�i�� ���J ��}��v��'Dt��V�-�b�d0F#�Eߥ�Q��nN���N��[��g�:�K!�)Tm:�G-�SS��t��c:�8���%�/�&c�
8�N�n���rRsS���m��\2 d.�k�����eM�<�6�8e@8�X2�����dh,�^���v�}1��Z]L��NE�f���F�n�q�q� s���&����+��-��:nذ��J�|��{Oߙķd�0x�V��hc��z��2���قN �$P`#e�MWL� �����b��`{�����}��}���WWD�\�rK��#\G�1��J�}�ABHI�X�j�r�ڏ����E��6�~M3��0��j��\��vF��B��D�8�hɞ�ô�>�G���&��x���Z�A/
���%vlt�ʆ�]�M>s������(�hK]�M�������s�"0��@�̌�i��p�>�X��9f$5���Q�� �b-���c��_��|Z>2y�Oh��#d4�r���4G�<?D4�%w�t�������D4o��=̭��D�b�*�����m	6H��L��
� (�us��ǦL@��#kp�pPp	��m���Ma�T�c��T3���G��%k�'_�r��P3@���g�w	U2�N%⫰NCA�B��4�@S���1�!4svX�"v2Q�by��������F��Mz�.��B (�l�V'� ��t&��
�2��������t\5R��E���n�↰�z�v;T��W�񢖂a ���L�;n ߯��,𶿄����Ȣ���_���	-�1�P����r'���¡�N�/O��] �Y�',IF�F���^(��Լ��P��x �[t�v��q B�<��f����Agԛ}�~/�L�Z;���΂C�}�H���Z*����&9�f������� ��&K�8&�ˬ����v6�D�tW�"��Y�X7�Vѳk��
ҳ7\�i9��{pj]�&��l���>�V��t�Ld_�*�"�A�Ryٝ���w�M���k��/����5	�7A�����i[0�}��cF�(Ta���Z`��5��eR��9d�"hq��%��m�ͥ�'�ro.?"L���w~�Yڿ ����T��9�嘲:"e�ػ��2�O�/,��*���W���o֒��'ˢ���e�T�s�`��ĭ�� ���#��p4O�8ݶ�bP hG��&|��o=d��b��d�H�V;	ש��}p��6�h����Cu���;��7n���-o� C�ґ�#�[v�M�:�G�#M��ݎf�ǘдӢ�����~��R�A���)�#E�VL�M=�Pke�j��~���Ė��-��dBtK�Ey���`�jC�{&S�N�)�k�%��:���0W��2i�κ�1J)te�ߧ�d����;�U�R�[�a�E�\�[��|l��S�ܟ�'J~�H?��s㱜��'[��`�v	j��	B-nC�F��O��p��} �g���?��^�����
ze�Z���Oˇ�ƻ���Y*⑏���`���a��~M���^���p�#���[)����&��,G��-{��!�X�VS䫵c����Z�Wqn3tdM�$&���_��7�`Y�2ZE=��_sCRB_���ePs�0f��;෿*u$�oB��霛�8��$7&�jɡi������N�SO��������a�� �NW�`T_l�.�Zd�9�BQ����T��%��K�j�B�&���6�x�)M�K.��6 mƛq���|�hh�A�x^��Mk�U}X�I��!�z��M��A�=H�N��X�u��>��HH.׉#ɠs�1�
0e{�:�\��ݟ�^:E��0��۞�&����Mk[]d�(nC_3?C��h3�S<����y�G�R����DH�ͤ₃[��{�aL���ķN���xؼ�#A��:e�TT����ҷY��N�K�Q���!����k�>r3U��C/��V��<֨��W~_-�3��fU:G���h��N���6�
U�����[q=f��<�Go����K8~��,�n��_�72dGU��[���1:_���WY-j�/푅UeR��T�m9׈	�x_0�x�"~�oNk����[<�� ������-�`5���x�����M	\'%?��N@$�Ϥ&t���<��Gƹ8`6��@C�f���X�S }O��I�/�;���$1�5-TK�dXZᾬlCyW2c���D��&Ω�a9z����2��(��R]�\ŋ닢Ƭ�8������g�궉�!�]O�XԊ�������7�
�l�M'�{;p�q�=�.r �.�։�'�JAq��6g����b]��t�QU�F�Y��c0�k�_�z	J8 r�5p��XT�9z�>u#Z"�@�/�̼u��f�w6Iv��ŊH__=�ݞ��c��������0��+�g�:AV�9Կ��a��%O@�X�Ot�U|H��vzG'p���\��Յr���*��+¾"������x�`F�ֽ��ϐ9�X�$���-A��
��|NO��v��
�%L1?O��ӌ�<�.�������w���CJ0]�u%/��1s3��h�51T���=1�ܞ�~�_ӣ�-�`�z���%�Nnn�	S�W&���1͡�R�^�w�[��}^˂�֍�b���Mi;gWۻ��%	mސ!>�0�#1ϓ����+��lǺ����<h(�ll?cBݨ#��R���8��L7]-����2psؘhs���q�b�5w���j9-U_�N7�
�	�a/ͻ�߷F����\IW�<J4�ICe�e�^��#.
�$3?��ի�c�+��5�bo��_��o��a��{$�	��Z�TC
���<H3F�]�C�D�:��3߮���ñ���)*��_��\�N�sx�`�)U��e���{*�>���wn��R���HM�p��B�`*��l��HH���}��8���l=��^���/:!�cQ���r�aQ�f�
�Z'v�|Q�ܾJk��F�u�zU^;U���en^�kt�u�bX^O��>t�������7���;A)+{'q�� ��ވ�}�b��;�,Y���I}
�l�/ʤR�5�L�}��������p����j���9���P�2�}�����dG���pv���3�@�w �yh�/�H�Ӎ��A�j�
Mjǹ���lm��Q?T%����8F��X�Eo��9c��r�Ko�o��xF)(��\+[�z���	O���ʁ�r\�0S�z��2 w��!P�~k@��ZEk������u���4P7��\�E�hމ�ԫ<�����=�
��:�� ÿ��LN��Y�6<�����	�W�B��hA���N�B�~}��C�_>�:�U��e�.��aig��L�t�����~��ԟ��ӂX��`�h���zְcڀ
7�����OC����
ƾ)=d�V9������0�E�is�6���?3	�%��Db>ps5FB@�)��ZJ���Z�j$T� �ھsTg䎷�`�m��:z%���R�j_�iJu��	¸��buh1gw���<l�5lA�d%�3 *>���~u#�� �L/_�b��l��CK�\u��SCq�O�|�K�:�8��y�!�c:����Bk?���o2��5��<z�'����fɑAa�ì.0}�ji�`�A7Ym��d�7D�W�]C}Q޶P7�m��Jߩ/1fJ��I��4�g7�+N@)��?�Jvb��upMHH���M��Jb�Ï�o��c���ֳ;6�R�{���K��՛ đ��� �`O��ټ����+�>R}%�q6�\Y~����ﴎ�P���L&��\�H>)�D�^U��-o�L�v�=�Y�Ԑ�;-u��9�f�ݎ$�/�Q���{,���D|��R�a�@I�Пr�v�i/��f@����R��P�ײפ�C��(f���?gG4�A���T5ڏ� ��ĮRb�+Y�6V��E�|.�V�	7+�}��[R�Њ/,1���`g��B����8 �ӂ?���"�.7��ȁ����M�e�0xgÒ�H!�ϕ<�k(�g]���9�
�UG|�ڐs �7�^�GS6!Qe���r6�$XF�L�)�i�5:q�*�OC�g���G���+cD��3=�u))�2c�(�8TU^j�ɀ�/�1��6�T-���Ha,XD������CFy�R�����)!�	:�Ք�����������Ja��9$a^��o��������W(�R[��
�����S����V �MrCm��U|$�B�m�'~�*���8��R�\��T$�wSg�]!���aS��	�EH���H��I��S9Z��ٗUf�����L�z`Ao[�JZ�fK�����`l�!qZ��t���|��*��������^!v���v)K%(L��ݞ3�&����m�T$2|T�Ú�|�/�TTN�6��g�M���~�W6�9��]3�w�\Xt�^���EǕ.��眑�	h����K')�q%����P�kF��g��T��8�5��gWq7gl�Ơ�9�$��C= 5�:o��[���9Y��#�!����` /��|]Ei)�K2�q<8o,N�����3�|F��)@,z�#���b:Jׅ��Cfq�� 3���V�[�d�{�����g�a23ڎf����eݏ�#���������?�/[]��9p�3 �JZ�>���ү��b�xݯ�7Zn}$���3�����N`k<���o�8O��s7��N�A^�̍|V�F�ĢT��6E�2��Ep�o�L'1�� �crmI�P���tS�ui�T����RP�%!n���Uk��B[���ಾ�́����)b@�B����Q�+��#�Uܨ����a��T�%k�'���sr�q��,� ɸ�"�P7��w����QI^n�l4��K��a=΅��j"��2�����0�( G$������X��+G�Mw�t7!cs�����ZO����9'���O���^v��:"�t��E
γԜ�Gr~��$�Ώυ�1�s�H{K�k�ܪ�$�Pg*���٠�2���G�	���H\#%�!����Cl&$�"T���������l��`^�qȢ�	�!�֐ce�i��`Ҷ�+�i���0��x*I˽{�Δ���#o�#�ͥ�i��F%��|��v���<����Ȁ�hl*��f��t0��Hzɜ��uF�O���O)��QsaU��i:�����F��8�؇ω1V@p�lC%m�ܟ^�z~(
4�G�j�0�,�qI*g;���r��� �V��9P�䖨�h��zOW�X�\,^�[�2��8�ǰR6�=N�©�ځI�b��W�l6���w/�u&�l��Nv�T`�׼��S�#�{�E�q0G���7]��0ďӃт|1J�Ox��U��&��аJ�� 0��W������}�3���E%���{��m���K���k��:�[�S�i;o�*�8��4�y�J���q�3x8$�W��#����EUI,�')��աwx�M\*�^��P���rm9�8d�z7__��O= �qx�L�D9K,7��*�`0�I�UoX�h!��`�9Q�MȳΈt��>M�aev�g��b;�ϞX��6�1���(z;y�t�J_������[��ۣ���oLLo�ېM��K?���@%��?� .v��H�Lw�L�
�Y�?����NQ�Lݏ�$�k`Fp��w�(�,ɵ��@��t��S�a�C�`�.@8d6[]DǨ����䌞	gB�zY/S�����[]��Ŝ]XUCF"v,+�{+SKX�V�C?V�u	�*�z՟�m�i`��%%ǻ47���m�<<+_���"@X��rø��/�L���E�<К�;u�����C�#z�׏!^���+B���]
T����-���A��x�2�
v!�8���R酮u(�g�.b�]�s�Tϝ���v=`fV�[<���N��b>�?��D���$TJw"<�y���zd������A0�H�ņ3�D{�o�W����z����Z��X�������H�>u<�T�"�a��6�R�c�b-
?�6&����X��B�I�N� �52�\8���*ٛ�����\�-��,�V7KV�n1@�8��X�
̜pG��a7�h�NV�����[:m�R�����#	"�� ��^J��-�Ƶ��P�oܧ��6)|iQ�M���E�2��ީ����8�_k��X����?%��l��uCD�'�iɣ��b��ܺ�c�2�6�h3��ڠ�D( S������ɬcJ�=�i���i��&��M9V��p��tT���<���4���D,tC�F��{��W�q�m݌vW7�N��C�/@)�M!Ҿ��{���!�"y��o���G.�����dn�I��F���ܰ�m�D��$�n��N_.MYt�<#���F�Dں�U�M6��/�?�A$f�K�ĸߏ4����th�폴x����c�.-�"[�`���,��w&���Py���7�Wi
��s'�k������0���@��Q vJ�\'
F��P(�֦���Qe.���{8��)�Y���	��̃�Y��2YpU)�1 ڲ q��*�n�Pٗ!<i,�/%?�׵�r�y�F�@B�r�$�s9{����rs�h=�ⶓYZ�>����	g��U�^�B��'�C ��T��q�hy(ɽ�6�E);ˤ|`	�����%W��M����L��<��&ڿ���^�U�����P?�^��hUR	)C&8p�G�5��D�C��)d�ŀ\�.k[��/em��UGS��	&bȨ��k|Ŋ��=���DN�%D�{���g[�Cvhe��jQf�{`B%��@�k���&+ �sp��4�.�y��#	�����$+;>�e��=����5��?���=��S�@^����a�M�b���`ӻjp[�F��JT/w�� os�6��K롊�P�%±ĳ�G+�"hI)�}ܡ�[-@�_��/�9��f*�Q�0`	�����z��V�g�_d&��ÿ)��=��|��S�K@F���i$�X!9���%�Zb +�9�i���G�u�
�nԷ�[gǛ�HM��0�J쳺����k��v踕a�u�G8kB�3���C�6n()$F�Na&�x��m��֨7���v��@hܦ��p������{<ˆ
l�\w��4��_� ?�%�V�.W-4Լ��3ֹ�=d��_�7��n��|K����y(Q���]�M�ɅQ�v��\}�y�3�4
S�)|&�X%���tW�r1J�w<���l�`�ՙ��������������\כ�u)Rq��1�*���w��P�n�v�a�9F�+������Q3q�[MuD%����j��	�Je����q�k��}��rB���&#q��~D��M�𬠸Ј����?�J	|2���y6��,В/B��yUN\C�h_W��o�OB�©J�#4�$�a!��A����ez8ϵ���x�0�ɲ.��>\TT�=�a���풑�]���Bݿ��r
��[�����>���j�4?�־�j2�����`I/����c��V]~t;��<�X�P<
�ˋ�@�����݃wV^r�p�O�zQ ����X,d�SMzx��&xFkQ���86��u~^g�ے�Z��de�0dR�l%G�n��}V��2�Gㄒ=t4�fBҁ#���15��.7UK�����?ɠV	�����'{u5���6�=�Hq;~q��{=s���~<��/��~+��w�eF~Z��ͺ�j��)��7��U�X���G~H��2���ɕ9��X����Q�v��� 3=�p\�2x�W8r؁v�4��el6�tZ���Sm�nSLğ�,��J�cߖD)̚ᳬd�(m�}~M]���?tߢ���?<n4Һ����M+g���e�Ksp4Ҍl�)&�wdf�cf�Ϳ����֪�+��I���e6W� d�r�u�A��,|��<���	�0�G��ͩH�p��:�ɑ�HjC�L԰�Yƌ��s��dc�Бa�r�Q����<,��0aC��Uu�X��@�w�׈�ڨ�#i�+p���ԃۛK2�Q�<j��eYۑIҜ������;P�mk���.q��Ҙ��L�h�st.����6ҽ��V�XbP��[��?BV���T?�bD��u�6�G�L%���lH�^����]^*�h��OM�IE��$z�Q�\�o�/-�̓e��f�N��#�j�g�Q;^��j>�0����"|5�V�%_�	=v���;�Idc��`�P�{v3z3��s=u�̒����P�����o��׌�Wy���qV������ǽ|�z���]�CS{�,[ Gs�X.��t�y�0� �w�0 �	�Rf������ι�$L��X��؉ND?�!��������/���%*�x�+5%o�ߑv�n:<O#=<]|�������l9�zN@��/+ˮ�D!����W�h�*���B�N{��������\]s����8!1�Xp�4[T��^r3�������1鷞'?���}�2�6I�'�)I���H�.N[�w�#=��d����$��&�um5���ܝ�l�L�y~0E����(��U��I���3��2���.�G����q�)p�J1ox@:�rEb��*�f��*��6P�ú�Nz�����#��t� 7�+}�2(���Ά��ݧ��S��[Aiﳞr���=o�z+���������nw��-�{�%���$qb�`�b7A��g@5x�੨�7��}`����@�ֶ�deUk��V��;`/�Bc^��.{�z�N>-%�>��N҉ ?�yDL����~o�?�������<��P*�_M����Jl5���(�&b<G�f�%N#M΋�����T�Ɠ���\��)���6Ѥ�ŲH�:]���F��^���K�PrI��k6T��H8S�ш>?���a��~QU��e��a�A ���w_l�Z�X�X+:YIa�w7��qYx��Y����U�{Ry��mj$	���ûv
+ 1P�k']��Q�+���̬��C��X��P�!���b����G��&ml@�nu�W>�,��O6ж�6K��J�L��"�����Z��T���Zo����7��1G�����JR�Rn��Ϭep�O���$ށ{��Q��lƘ�����`�k���~S����dŃ�Ÿ��#��B��K k����՝�WB����}�p���]��M�6Q�KU�t׉h�ke-HC�+��)�m���|�x���,m�ɬ�%�rũ���?���qQ��z��M�Ǆ�����^m�'\\&��2��,�i��T#�t.0���?��_��[��ӆZ=�lyC�n9OuiM���[��>�D���a��f�ʡ���6�lo�j��XآK���n˜:��A����)�pN��R�j���,m�MU�ȯm�/�,�CvU����/w�l���L�X������v�U���M"��~���b0�|��d���r��Nr�2���k�}Y��	���m����dh�b��j!��5l�.�' �i��F��mݺ�1i�!$A�.q���M��m
9�;����	kw�I�c�BZw������DP����%Gq�b�&3G�,���Hs�@�4߿ )ŏD.)�A� !$I~TnB���L�GFԷ�uaު�W�@qlh�e�i�h����a�N��β�M?G����<E&������!��av6�Q�7|Q:)�:�'��3���%�W��g���i�9X���N3�혚��+�aX�?o��c�ݭ�bN�R�[���ͷ4�{��{���v��Z���b�-���盄M\*��e���ۯN��D҇�0����s:�/%�(7��Qݙ�<KC�ђ���nw!�����RuH_�����f��t�˵�O[�e�?�U!F�o��aPb�ky!e%���0���������on�LH��q���<��œ���`��ȅ�l0w�)������UB��L<�j҄��q�/��\4��N3�\��#�-f}��g��'�(��Rl4�R���q�~�uz?��0�U�5S�;������.O��"j�~��L�0xT��E�~~�i����O�,��w���lYӈC��9}e3�X&[Tu��{vvKl�_�3�Ƒ���R�0Gc*�B���8�P͝}oi��V-���_m��3���b�
�S�^3�-͐��C�\�:Di��p5�h/�I����C����ܐ`٦����]PT@u���0�K��)m�M�%z����|��ۢѹ?�G��_����L��L_GP���Ʒ*$���Z�'3M��!D��#k�]_F��%Ǖ���U����ս���)�1V�_:�QEGru5��h����y*�����PHs,�a��R��(ƌ��1�
���Yz��Z���a�� ސ9���d-��.����	�qr �Q+[���8�Ԗ^�������zepm���C�r|iS?��ғ�Q�^V�d����]&jM��"�j�oY؀�'�S�5u�����)=oq�Ҋ>��On�������GKxд��y*j`���G`!��|_�.��!�b/�m��BH�v��6p8FL��W�o���3_M�޻�:Z*��d�W�k<x!4$�57,]<��
�oi"Lr��+O(�Έ
���M��ʲ�2��Rc�p3`��[�ǫ:r�l'���t�����RZ�j�����ȹ�֤��`鞊��~*�"��s�wѾ|6u�x|%i7��e��'�V
$s��B"D)?��*!5����IG�v	B���E���>��:왯������R�΁i�G~SN(I�t�.n�C�LjWЖn�����e��*��sѕk�~<&���P�l��f��U��SM�Ku'|���;_���zG ��Mɖ6�iSd���Q}dL�}[UY��L?�􀢹$���͌��%�ynU�5f� � f��Feg��g�I�c)�ja�7͂aUU8�R��'e�D�r�F�n>�C��A�KX����=�2�[߄�Rכ�����1v��p<���/��h�����p�����բ��ɢ~W�I�O����f�`ә[c���pT��E��8�I��:�1�A-�rw�*-�����j�"�.��"�N����ͨ��cS�h��9<�Қ���5�j{W̑�d��d�N�p48�P,C�����'f���֝�"�9(�����R嬆��OG�2N!x�HQH4������zG[����2���7~x O3��H����x!������|� �-p���s{��(�RJ�^��3It���o���R�@s5����6����c9�lk"����^������ݡ�6��ŕq�I��d^iB�"ʼ�X0Hý���S�`>��B�wo����f=wG�%�!�_����k䨜�;�y��r���b��u�nF�s#����-i���@{a,6#������v�h�&��``�V��h��c,�#��giX�z�w�2/�/���}�f�WL��:JmMqw��Yv�ǑFmҔc�~�^ҿ�>+�VM�(Y֬�:'����e#��wK3��EN��f��ܾ�%qB.B:�BOn�q��O2����]�r�?��k��9��E�C�r#�J� �˅�Վb1��mN!��N���"U�!,���y������������$�U}���X��'�Jy��N���|0F�A�oR������)���x�N��͠��ʴ��K,�%�ۖ~&T�7l���k�2�C0��U1����R΍��<��a��f���z��!���,?��K��F��ޭ涎~�%�7�/.�j�=D���#0�)�E8��V#-ZtF�];���=�LE3�/8c���?F-_�+�PZ$t[FQ�6�'`���}����E�9&�R�w�@�`��
!rt��ٝu�\ҫ�ݰ!~w�x_���"�s�H��FZnGh3P��'��ڀ�q3BܘSF[-<���gK����2�Z�$�ܔn���Cwcj/�_�_Q��0�"J�lP�&���
�wF��"�x30៫�.�8���@m.�,�ˁR����l�ҫ�����?@�m8:��3�fY5Z�:�>�
@�X����� Y���}�n�HDb�uem����g��42��m�����h8��)/���DY���x���*EQӻ�r m@��/*��6XSeX�i�D-!�'�|ɣ���n�Ѳ�[�9p��Bb���ǨJTl�L^�zl�&ԁ�Jz�;�S������z��C}��f+-zLfqF�_Έ�;�.Czf��}~^�������6C���2��E��7/���	H�x$�g��8��5�-8o���A��}��R�ܫ7mp<;9N$�c"��t���C��;�r����N �d�Xa���k�S��,�s�zM">�t����a�J��s�H1S��n0z����K�Qym�F�a��`����ŊNB�?wI�3��tף Zl��~f�g��e���'�H��WH�oX��@�=	=1��Cl��S�i��E�>Tׅ��ͱ�R�T�g��[�1�P�� ��M٢k�P��tK�K��p/Bx�v�C{~�4��~�j8�-jS�_
WbEEb6⎺���I&��I���'v�e��\�s�c�Z@�Q�����j��ܔħ�@��3`�4S�8V)HwG�㥊$�6l��W�A����|`�*a.���	������+�n���МEB��t����lGQ�1=�UU}S���&v&���>*7�^$M��N�_�c�qޗO�������[���}Jshj���S��K�c>q�]���v�:�!A��s��7'�U�z���=�� ���Y k*�/����o��C�'K3f_� $&d��~�b,!��%Ź�ֆ�elJ��t+���Q+M݊��ip��M=�)��3 ^�g�Y�v�a�6Kh�y��I����"_Wu�%�E�H4Χu]�)�-'�)f�'�b��o�&K��J3����*&ړ�Ms����ǉ�ߠ����a$�鯘��Axc����S�6Z��0<�
+�zc`4��U�=��LS�Y�a��@2 ��bޑ��4�����@���o�j���{\���AN>�	]�/���md��b��号�y?1�����1���	o�L���M�{�ʕG����o��({�g7��"	�Z�v4����ߙd{;��-p��[,�bR�[b���F��ɕ^e=�?d�1��v�z`z����g`wT�d՜�h�2X�o�h;ꉳr��Ĭx$��4�<v��N��͔�qr��)C�Wǹ!>}M�/IY]�!ø�f&�m4(9��aA&� 
�y�M�H:t+�o.�y�{к��V����N�ڿ<x=�;�1g�s�$��Z�Q����h4o~��P�a͇j�^	��\�oDOB�=�Ϥ��ٽ˥�v���4�qi�?hC�m�?vA�X�J���Ws�3�����X�s}��&C�@A6Q�A�iU9���_��� �M��^�բ�����j��r��+�br��8N�7N�_/�T�w�s@9"{��s�k��~XwY#�"��c �?�Ζ�!���E���1�(*���W�V�cǺ.�?�G�����2;h�������\��.m��.a#�J��Q�|lS��i�6d_�_S�QլAfڗ����B�ah~bu�mP!VP�/u�mO�~JF��_ �nI^
�XdC�]�=�I�w [���4W�N��������)�v���].�a��K�D�BK�\,6J��hx%����J �C{�׾�͉S��IZ�Q&W ��,;*T�85��:������Sh�|�ʵ�`:E�g���-A2'(�C���b��ت����c(��9�Q��+�ٯb��,]�������j]q����0��B�!l%=1�`��Qj[{t;/�[~�<�J�r%�=W��������p�'��r�N�9>�
F���yO�1gm��C���/5��t	#-o�!ꇻ�($b������5��?����5�϶�B�p�#���ŏ���
���j��}�x�w��L=�w��9��)���̥�\�im'��2�y�܁s�����F��n��?1�;V�kp�=�]�S��M�t�"�;N�z hg	 ���e,)�!�J.�M��`~��.ʈe�laH����MW=������`��ta;-��(�j?�!'UI�h���p�8��GX�nY�?�����]��TI�u�"�+��_�pڎ�=O�J��.�4.���8�z��#J�P���r �Q���>@���b�͟e+��e��P1�~WHi>3���z�=�˦ş�Kv�4��B��#�ҕ?�)���P٣�!!�]h��W�`��9^jn�af�s�|e�)�W]�d�KM�8�*�m���A@摔��`���h��q_�s�0�r��������;�j��7u?L>�|����2۞�	�ѹ�������U�.�eT;#�[q5��و����� ���-���-�=a)�ާfD�����v̩����8����m��@ǰ4a@����?y��5�?���QR.�+��+�i4�sP4+G.W1���14j����p����*2�3*VGՐ���C�}�O�A��r�@�'|xn0�<g�h2N w�v����_� �m5`�pbM}�@)�~(�) ,#C���Z���
�e�4g�K��wv�Idm��g2��@���4\
�~�`Cd2��`��������م8� �I*������nC��a�@�����#�F��,�K�����w�xD.�m&�����
�T�)��H)���H]Yʻ��ԝX[�� �*R�?X���5��ł^�
q�~)z�]25b�1b�*PB�Q���K�ρ���3�ɲ4#ȆY�JUԳW�������u�����GL�G�����L�îrw��Y�f�x�U<�~�}VL`���\�Y
�<�v�:Z6ƌYd9-u�aӸ��s�ǋ(_��q�*�$��>�v)����E3J�*6��ڣP/J��������@? �a���}W�R3��6��W"�Ѭ%\�
�g���-���j�TȨ���S��2}�a�t�:��EPm(�<-f�X�mL4�{��!�+3h'��X��p��Fq���{��@�`�b��� �5$
_��^�g��Vd����~����ߗ�K���Kߨ��UЭ;�\��h��3V�b�Np�'�[Kf}��F������F�����=�������P�au6�p��A�����e�@Ӯ�ūO:,��ȅ'�����s�)��G����{�ܦ?�ÖW��a@�/�g�)�ۍ��`Mmj���ۏ�{�z��E2SłD�;�- ]�ȼH������B���b��Z�����0�,H!���E�T
��|@8��f�"��|8Sl-��xhkT���^�㘠���Է��7Cw�4[�����D�����wsq$v�c���d
�L^J��v�R6.���+/��X
:�g�Vv�S=��Fe�+h�D{f@ᅝr�û���O���(ٲ|���e D^ثqvF�����ab�x�����C����[�i�|p<���K��m��sg�ԟF�/W�|Y�}�Ф;e��P�y����#�EDY����ƈ�#2y���"~�x�xn@zn��b���`,�0���Y7 ���k�?�����XG%����Z�kE���ڤ�#D�'�j�^�=���e~}_y��~E'!B��{�i�Z�B>�$pO�e�ʊ)<` ��5����\_G��I����w.Eh�g�L<*(��!�H�G�F��捱{ܝ���5@�6(�}�4ѻo��]ZY>��f�b��3p��Nn{�/X��L����Y��&	u����V���^%��
���B�p�[����]����F3m������u��ra�4^3����%\�Q�9��Hӹm�?�О�V�IA�&&B�GT�����\����f.۷<���R
���CDg�0w/����0���	�Q����BĴ�Lą6����w��$�C;��&���C�;�2���n甎0W�~ԗ�+�Q�f
�Ƿv��%�9�������Z	Mٰ6��d�S�Ǒ��3rW{E��a���b�I��@��Tj)��=����C�M�� �����IY �K�Ȉ�X1��7T#�\�N2N|$���{Rc���ke�X�H.��E�.o�/4�ɢ`=q�|�SH��ϛ��$n�+Z���>RF��K���]W��ї����$��UKm��.JħO�l"8��C�y].�2K�����h�S�y���x$����+�!�{W�j�	PU{^qQ[	?��Ӗ�Ǌe����L�o�I�]�l�g݀v!N�P+�� x�!�P�����ߚ$�#8P�}��E�wk�vv��'��� (��M�cTsCky:��*��'��?��L�=�����/BЮ�J�E�|���і�:��?e
��A�#O�@��0�i2Цi��x�X,1m)B�e����/��rs����u8W�ct���|9��9��@̞�Xԕ���a�d'u��"����}bh>�ݶ`�v��/�Xڶ5�!af���|0�&����htר.�DumlA�H�
�;�E�H}6�$p;���w���u�f��I���SK�D
��\���/{-�"c�5Yװ#���9�`@׻�ff�9�J�WL.��*�XQ�����S�s%����#���D��9�S4��@���xgff,��e�����lzˑ?�[���w���u�أ�z�
g��7�Q�)�$�hy���$��U���.��Ϥ�e	x��.qN�l��q��r�Vf�Y��k.E�U1�d,�,����\���:X:�^K�&��	�Ƃ�"�̾lٹ�� �]D��q���5g���ݤ�����z�T�
�_��	�!x�6��0$���t�q��,��������[+;ˍ��A�8X}��+�����S�����t��du(�e��z�RC%���G���t�Ij��&�t�;,�m�K���1�;��'%�w�ʼ�&B�?֋$;e�]H�/�RM߮Lj�j\]������}� �'IZ�-Ep�����Jև°*ᯑ�h}�~?O�H��g�V�lԹ�� �{�6O����Ha���8�A!�I}X�C��LU\9�H����$�1��l���(�lnz���s��/>.�	q#$l�t:i6�d�&J��*V�__2H��=#S�(������V����5�|�I��"���$��hE�|�����gAgY�W�fO ={:uW��D�Pu�����*}���M�cD���:kъ]N��?��z�L�4�Y"]tiǆe�+1�B'��rx��=o��!]�IT_�j���{eç)��A�������+}g�XRx��Q�H~��$����2{�,G���a�2"����$�a�lT�Oh>h��\�DN�GUqs��*��t��{⚟��L���yFV��y�͍�#M8�+�ۻʂ����׵�h6(Ygո��4VK��M�V"偗�N��N��A�T��[��R.�(�3�L�u�L�ڣ����]qIRV��A�U)����6�J�G��y��O`^��}E���IVi�	o4�7s#���kU����`�G$s��(Ņ���K�J�[�q3���I�S-�'�5:�@0г&��l۴O�GJ�s�>�"`��g��Jy�����n�g��LE�T�
��\��A���v#��5k!��I�A���G8������;�����#	bT.4��DUv���?��0tgq��Q����TX�(��k1�6]�e�v���z�nN��:N�����Ot�eS�^-�m�M�g�K��w�LR���������2�Q�}�����Ml��ِfIG��o�Y�u7i��' �(�!~!����@�H�r0r7;�S��m@�D&�)��Y� �|�z�sWv��Ÿ��÷ŋ�D���Z�����ЭT>L�����Q�ow䉠yC�;%u5���}�*ӎ� }	/�#��(C!�U���G��Oi�5���S`\c8 �s7�FrT���q��۩I��ȢvD<������a@���g��^/�j�;��t_n`ݭu����� ��d`�{@���+|]�R��'S2ˮ��ߣGvbrm~k1N���T��d�����BIN:�	�tgq��E�u4�����������F�~1�
3���w�a� ~�}��wJo?�PҔ���F[��*��[�IU,i�Rq򀸍�6�� B�`�~�֥V� a���j��`&�Tq ew����w���h���պ�)]��0^� �2�T
X�-�*%�'��`��Z��!=�9����qq�����E��)bP^����o�/���:��5Py�M{��h��M%���P�n��F�.p�|�6:���ӭ"��@ ��hX7�2��YzV/��r�Ƴ�[��F��j���J�B���苇���ŝ:'�m��f<����Ֆ�d[\p�5�ޡB�ߤ�u�~Xb������R~���S����w<:��9�}�uN�-0  +��O��w���ctue�Q�*�ܖ֤�9n��$�%M�4l���d�a.\o	���&����^�arBz��(C@l����#).[�)�p!�,Y�3N,,����t�(ia`��J��Afi��8T�N��_
g���'���),��P%�>C˺��K�v�`)fk���G�r�g��G���m#T���&W�GY%�V�Vd���i���E�ǹ��x+m�ڪ�j�+�J^3;�}C�*S)�d�oN���J:��v�j���JI��i���8=�����08�+������)C_g���C�Q}�00���Ի�3i��{�_��Hk��;�2��O�=����|�P��6�`��`X��	C4>�j'��O�d�<����`)͋I8|�Vꞁ{ٽ�EP`C�D�Zzf�E͍�`�4�#��E݊��%Ȓ%�N) �]cK�!�=Zn�B_FR��F
�s 	/�2���$${�)/�׋����Kب�+�j�Y
��
����3�
;jj?E��ӷ��OLT��JS�Nk��]&\�fB��,.�B%��co��U�k�c8|H��H�\���u-/��j������<�g�Ub������fݜ�H>�V:�z�����J�鷩I�I}�d"��F��;֫���	i��SS'C�7��E.6�P�SGP��t��ɐke�H���8��`�X^�*9)C)�4��ǈ�ov�,��U���5�+3wy��pϤ��ߵQO奘�j4�M�MV��O	����}�>#hF��k�WD-�%��'pFm�s����F��˂��0��,��ss%I$�cZ��H[uk�9��p�dv��m?0����κ�GۣR��٣�����7�>D=4�̀0z��v���V=*Yx����Z�?�ɯ����vT���ge� �q����^�(�������y^f�4�ʏ��8>;�h�v,Bd�Pn pZ5����,�EF�`z��`��]���fcg�w��`';#m��g����Y`@a��*  _����m=���ѯ��؈9�eՆ��)�^V��� '�$�|T\��1/rw� ����Bx�ӱ@����I�K����O�J���_`B%p���ҕ� ��@r�zT%�0�{@l��x1t(��uO߽0��f����V��Ħ��i�,��|_g�~�&y:��>?�G��ԴO|+���g#lA��(���K4��4�˖�1�9�~l U���~�b����H:k�|p=����8y9�g���Ԑ�Ti�3Uه��9҉�<�?�2��M����'y�`I�%����i��Vt��9I v����Iɇ�s�C��o�w�TaL�׬�f�� j�6�4�q��[�'�V�=�+����"��)"���k��8�~)U��:�R4S�g	�	lm����� �|�B�7�����2�i�FB�a��E8�}�M���;e&��rm��Sk*.7NA����p�k����O�VӇ-~;���n����/�0"Em�vNwI���/��1�lb�zZ輙����TaZZ3k��_�{/IQ�V�Ѐ��]q��,��Q�9�3M��hx:���#P �Ya��Q�������Z[R��>��G��.����8��;Xc��B�Il�F�l���;ԍ(�DK�)*�^�lL2�ɜ�3��"� �t��1{�h�Â�Q����x�[s�CO�����p�Z��"x�QJ�gpv F6]�H:�m.�hP��Z9��j�$k�J�u� qfI���_�R�_~'q�PX�t�кԪ��L�ec@e��r�p[�����|[�pٷ��-��E�dk��s����}���پrHҊ�E�������l�J �B=n�$��	�n	���%���"���l`�|I�,��:�J'%�=�����HcO<��#P���Խ�T?K�/-��S���t�w���nK#p#���P��Ys�"a�� �t�P�+4]ui��@�\�vu�d�?�A��� `���o"<��4�-��ǥr��D�����
Y&�f�T}[��.UD�H\*��|c/8�p̈^��'~Æ-�۲�˸���:sJ�(ua����j?�(��S�!���)���~{���$���H��U���}�W���*&&AQUfwڷQ?u-5�����a�ρ����ZL݉�����k�$Mi��Ӥ⹎}wQr�x�NWݵ�h�V���`MB�H�6�\y�2N0���;'M�Y��ó��V���vb��'iǄ��EwV�L@xy���[ZY�OC{��]�Y����l��Z
�1���*�F�eR�|hH�\��P��@�=J8(�����[˽Qs����_��S��1�J��*�6aS,�������$;(ؕ�i����Z�<��Oi
0M~�
�փ����cF�U%D�������C��d
G�����;I��0Gߺ�0��?{���x�fn��6�82��h����q"�V��&�)j�H5�R��8-�3�/ޓ"�,h;�:Ӵ�����P�{���x�FV�AXm����c�8��^	�\�����2&�G�M���T_�D�Ω�Մp���:�a�Z�2k[��QQ��d��^��*+�U �^����������^�]��$�|��J�ܢ]�su1����՛Po�ˢ���O�{�=�q�$��e�)sv�i����f�7��t�J'��/PV;�V&��0��>�jC�1����@���9/�u�.⫾�:��3�E4�J\�*T���0e�a�e����^���(?Ukg�|�}��]$�ig����"آtN�ַ�f��Z�u��^{W�+�'~6C8_�,�u�ut��W둺=Ϥ.��+�BQ�7��bMj��[���?:��G���C��L�E^�R����G�]�S�륹�_������J�{��;��s���W>h��7�4�/l ��2b���/�W�r�N�.�=