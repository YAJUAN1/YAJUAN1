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
  static void plotVecTangl( int col[], int coln[], int o��&�	U4A�
"(�S9r����_
�9��X�[����3%���*�,*�"���䇃�jĲP2uL�2�(b��� p�BBg�k�T���
��,�>F�+4���
�8�rߖ�T"aS��R��'e�1�bs!�z='9"�e?���wU�����v�w�<��{oId(=�z�0ڔ�P�ą�A�����2�ͥ.Eͭ^�}�l�i
8�N�n���rRsS���m��\2 d.�k�����eM�<�6
���%vlt�ʆ�]�M>s������(�hK]�M�������s�"0��@�̌�i��p�>�X��9f
� (�us��ǦL@��#kp�pPp	��m���Ma�T�c��T3���G��%k�'_�r��P3@���g�w	U2�N%⫰NCA�B��4�@S���1�!4svX�"v2Q�by��������F��Mz�.��B (�l�V'� ��t&��
�2��������t\5R��E���n�↰�z�v;T��W�񢖂a ���L�;n ߯��,𶿄����Ȣ���_���	-�1�P����r'���¡�N�/O��] �Y�',IF�F���^(��Լ��P��x �[t�v��q B�<��f����Agԛ}�~/�L�Z;���΂C�}�H���Z*������&9�f������� ��&K�8&�ˬ����v6�D�tW�"��Y�X7�Vѳk��
ҳ7\�i9��{pj]�&��l���>�V��t�Ld_�*�"�A�Ryٝ���w�M���k��/����5	�7A�����i[0�}��cF�(Ta���Z`��5��eR��9d�"hq��%��m�ͥ�'�ro.?"L���w~�Yڿ ����T��9�嘲:"e�ػ��2�O�/,�
ze�Z���Oˇ�ƻ���Y*⑏���`���a��~M���^���p�#���[)����&��,G��-{��!�X�VS䫵c����Z�Wqn3tdM�$&���_��7�`Y�2ZE=��_sCRB_���ePs�0f��;෿*u$�oB��霛�8��$7&�jɡi������N�SO��������a
0e{�:�\��ݟ�^:E��0��۞�&����Mk[]d�(nC_3?C��h3�S<����y�G�R����DH�ͤ₃[��{�aL�
U�����[q=f��<�Go����K8~��,�n��_�72dGU��[���1:_���WY-j�/푅UeR��T�m9׈	�x_0�x�"~�oNk����[<�� ������-�`5���x�����M	\'%?��N@$�Ϥ&t���<��Gƹ8`6��@C�f���X�S }O��I�/�;���$1�5-TK�dXZᾬlCyW2c���D��&Ω�a9z����2��(��R]�\ŋ닢Ƭ�8������g�궉�!�]O�XԊ�������7�
�l�M'�{;p�q�=�.r �.�։�'�JAq��6g����b]��t�QU�F�Y��c0�k�_�z	J8 r�5p��XT�9z�>u#Z"�@�/
��|NO��v��
�%
�	�a/ͻ�߷F����\IW�<J4�ICe�e�^��#.
�$3?��ի�c�+��5�bo��_��o��a��{$�	��Z�TC
���<H3F�]�C�D�:��3߮���ñ���)*��_��\�N�sx�`�)U��e���{*�>���wn��R���HM�p��B�`*��l��HH���}��8���l=��^���/:!�cQ���r�aQ�f�
�Z'v�|Q�ܾJk��F�u�zU^;U���en^�kt�u�bX^O��>t�������7���;A)+{'q�� ��ވ�}�b��;�,Y���I}
�l�/ʤR�5�L�}��������p����j���9���P�2�}�����dG���pv���3�@�w �yh�/�H�Ӎ��A�j�
Mjǹ���lm��Q?T%����8F��X�Eo
��:�� ÿ��LN��Y�6<�����	�W�B��hA���N�B�~
7�����OC����
ƾ)=d�V9������0�E�is�6���?3	�%��Db>ps5FB@�)��ZJ���Z�j$T� �ھsTg䎷�`�
�UG|
�����S����V �MrCm��U|$�B�m�'~�*���8��R�\��T$�wSg�]!���aS��	�EH���H��I��S9Z��ٗUf�����L�z`Ao[�JZ�fK�����`l�!qZ��t���|��*��������^!v���v)K%(L��ݞ3�&����m�T$2|T�Ú�|�/�TTN�6��g�M���~�W6�9��]3�w�\Xt�^���EǕ.��眑�	h����K')�q%����P�kF��g��T��8�5��gWq7gl�Ơ�9�$��C= 5�:o��[���9Y��#�!����` /��|]Ei)�K2�q<8o,N�����3�|F��)@,z�#���b:Jׅ��Cfq�� 3���V�[�d�{�����g�a23ڎf����eݏ�#���������?�/[]��9p�3 �JZ�>���ү��b�xݯ�7Zn}$���3�����N`
γԜ�Gr~��$�Ώυ�1�s�H{K�k�ܪ�$�Pg*���٠�2���G�	���H\#%�!����Cl&$�"T���������l��`^�qȢ�	�!�֐ce�i��`Ҷ�+�i���0��x*I˽{�Δ���#o�#�ͥ�i��F%��|��v���<����Ȁ�hl*��f��t0��Hzɜ��uF�O���O)��QsaU��i:�����F��8�؇ω1V@p�lC%m�ܟ^�z~(
4�G�j�0�,�qI*g;���r��� �V��9P�䖨�h��zOW�X�\,^�[�2��8�ǰR6�=N�©�ځI�b��W�l6���w/�u&�l��Nv�T`�׼��S�#�{�E�q0G���7]��0ďӃт|1J�Ox��U��
�Y�?����NQ�Lݏ�$�k`Fp��w�(�,ɵ��@��t��S�a�C�`�.@8d6[]DǨ����䌞	gB�zY/S�����[]��Ŝ]XUCF"v,+�{+SKX�V�C?V�u	�*�z՟�m�i`��%%ǻ47���m�<<+_���"@X��rø��/�L���E�<К�;u�����C�#z�׏!^���+B���]
T����-���A��x�2�
v!�8���R酮u(�g�.b�]�s�Tϝ���v=`fV�[<���N��b>�?��D���$TJw"<�y���
?�6&����X��B�I�N� �52�\8���*ٛ�����\�-��,�V7KV�n1@�8��X�
̜pG��
��s'�k������0���@��Q vJ�\'
F��P(�֦�
�nԷ�[gǛ�HM��0�J
l�\w��4��_� ?�%�V�.W-4Լ��3ֹ�=d��_�7��n��|K����y(Q���]�M�ɅQ�v��
S�)|&�X%���tW�r1J�w<���l�`�ՙ��������������\כ�u)Rq��1�*���w��P�n�v�a�9F�+������Q3q�[MuD%����j��	�Je����q�
��[�����>���j�4?�־�j2�����`I/����c��V]~t;��<�X�P<
�ˋ�@�����݃wV^r�p�O�zQ ����X,d�SMzx��&xFkQ���86��u~^g�ے�Z��de�0dR�l%G�n��}V��2�Gㄒ=t4�fBҁ#���15��.7UK�����?ɠV	�����'{u5���6�=�Hq;~q��{=s���~<��/��~+��w�eF~Z��ͺ�j
+ 1P�k']��Q
9�;����	kw�I�c�BZw�������DP����%Gq�b�&3G�,���Hs�@�4߿ )ŏD.)�A� !$I~TnB���L�GFԷ�uaު�W�@qlh�e�i�h����a�N��β�M?G����<E&������!��av6�Q�7|Q:)�:�'��3���%�W��g���i�9X���N3�혚��+�aX�?o��c�ݭ�bN�R�[���ͷ4�{��{���v��Z���b�-���盄M\*��e���ۯN��D҇�0����s:�/%�(7��Qݙ�<KC�ђ���nw!�����RuH_�����f��t�˵�O[�e�?�U!F�o��aPb�ky!e%���0���������on�LH��q���<��œ���`��ȅ�l0w�)������UB��L
�S�^3�-͐��C�\�:Di��p5�h/�I����C����ܐ`٦����]PT@u���0�K��)m�M�%z����
���Yz��Z���a�� ސ9���d-��.����	�qr �Q+[���8�Ԗ^�������zepm���C�r|iS?��ғ�Q�^V�d����]&jM��"�j�oY؀�'�S�5u�����)=oq�Ҋ>��On�������GKxд��y*j`���G`!��|_�.��!�b/�m��BH�v��6p8FL��W�o���3_M�޻�:Z*��d�W
�oi"Lr��+O(�Έ
���M��ʲ�2��Rc�p3`��[�ǫ:r�l'���t�����RZ�j�����ȹ�֤��`鞊��~*�"��s�wѾ|6u�x|%i7��e��'�V
$s��B"D)?��*!5����IG�v	B���E���>��:왯������R�΁i�G~SN(I�t�.n�C�LjWЖn�����e��*��sѕk�~<&���P�l��f��U��SM�Ku'|���;_���zG ��Mɖ6�iSd���Q}dL�}[UY��L?�􀢹$���͌��%�ynU�5f� � f��Feg��g�I�c)�ja�7͂aUU8�R��'e�D�r�F�n>�C��A�KX����=�2�[߄�Rכ�����1v��p<���/��h�����p�����բ��ɢ~W�I�O�����f�`ә[c���pT��E��8�I��:�1�A-�rw�*-�����j�"�.��"�N����ͨ��cS�h��9<�Қ���5�j{W̑�d��d�N�p48�P,C�����'f���֝�"�9(�����R嬆��OG�2N!x�HQH4������zG[����2���7~x O3��H����x!������|� �-p���s{��(�RJ�^��3It���o���R�@s5����6����c9�lk"����^������ݡ�6��ŕq�I��d^iB�"ʼ�X0Hý���S�`>
!rt��ٝu�\ҫ�ݰ!~w�x
�wF��"�x30៫�.�8���@m.�,�
@�X����� Y���}�n�HDb�uem���
WbEEb6⎺���I&��I���'v�e��\�s�c�Z@�Q�����j��ܔħ�@��3`�4S�8V)HwG�㥊$�6l��W�A����|`�*a.���	������+�n���МEB��t����lGQ�1=�UU}S���&v&���>*7�^$M��N�_�c�qޗO�������[���}Jshj���S��K�c>q�]���v�:�!A��s��7'�U�z���=�� ���Y k*�/����o��C�'K3f_� $&d��~�b,!��%Ź�ֆ�elJ��t+���Q+M݊��ip��
+�zc`4��U�=��LS�Y�a��@2 ��bޑ��4�����@���o�j���{\���AN>�	]�/���md��b��号�y?1�����1���	o�L���M�{�ʕG����
�y�M�H:t+�o.�y�{к��V����N�ڿ<x=�;�1g�s�$��Z�Q����h4o~�
�XdC�]�=�I�w [���4W�N��������)�v���].�a��K�D�BK�\,6J��hx%����J �C{�׾�͉S��IZ�Q&W ��,;*T�85��:������Sh�|�ʵ�`:E�g���-A2'(�C���b��ت����c(��9�Q��+�ٯb��,]�������j]q����0��B�!l%=1�`��Qj[{t;/�[~�<�J�r%�=W��������p�'��r�N�9>�
F���yO�1gm��C���/5��t	#-o�!ꇻ�($b������5��?����5�϶�B�p�#���ŏ���
���j��}�x�w��L=�w��9��)���̥�\�im'��2�y�܁s�����F��n��?1�;V�kp�=�]�S��M�t�"�;N�z hg	 ���e,)�!�J.�M��`~��.ʈe�laH����MW=������`��ta;-��(�j?�!'UI�h���p�8��GX�nY�?�����]��TI�u�"�+��_�pڎ�=O�J��.�4.���8�z��#J�P���r �Q���>@���b�͟e+��e��P1�~WHi>3�
�e�4g�K��wv�Idm��g2��@���4\
�~�`Cd2��`��������م8� �I*������nC��a�@�����#�F�
�T�)��H)���H]Yʻ��ԝX[�� �*R�?X���5��ł^�
q�~)z�]25b�1b�*PB�Q���K�ρ���3�ɲ4#ȆY�JUԳW�������u�����GL�G�����L�îrw��Y�f�x�U<�~�}VL`���\�Y
�<�v�:Z6ƌYd9-u�aӸ��s�ǋ(_��q�*�$�
�g���-��
_��^�g��Vd����~����ߗ�K���Kߨ��UЭ;�\��h��3V�b�Np�'�[Kf}��F������F�����=�������P�au6�p��A�����e�@Ӯ�ūO:,��ȅ'�����s�)��G����{�ܦ?�ÖW��a@�/�g�)�ۍ��`Mmj���ۏ�{�z��E2SłD�;�- ]�ȼH������B���b��Z�����0�,H!���E�T
��|@8��f�"��|8Sl-��xhkT���^�㘠���Է��7Cw�4[�����D�����wsq$v�c���d
�L^J��v�R6.���+/��X
:�g�Vv�S=��Fe�+h�D{f@ᅝr�û���O���(ٲ|���e D^ثqvF�����ab�x�����C����[�i�|p<���K��m��sg�ԟF�/W�|Y�}�Ф;e��P�y����#�EDY����ƈ�#2y���"
���B�p�[����]����F3m������u��ra�4^3����%\�Q�9��Hӹm�?�О�V�IA�&&B�GT�����\����f.۷<���R
���CDg�0w/����0���	�Q����BĴ
�Ƿv��%�9�������Z	Mٰ6��d
��A�#O�@��0�i
�;�E�
��\���/
g��7�Q�)�$�hy���$��U���.��Ϥ�e	x��.qN�l��q��r�Vf�Y��k.E�U1�d,�,����\���:X:�^K�&��	�Ƃ�"�̾lٹ�� �]D��q���5g���ݤ�����z�T�
�_��	�!x�6��0$���t�q��,��������[+;ˍ
��\��A���v#��5k!��I�A���G8��������;�����#	bT.4��DUv���?��0tgq��Q����TX�(��k1�6]
3���w�a� ~�}��wJo?�PҔ���F[��*��[�IU,i�Rq򀸍�6�� B�`�~�֥V� a���j��`&�Tq ew����w���h���պ�)]��0^� �2�T
X�-�*%�'��`��Z��!=�9����qq�����E��)bP^����o�/���:��5Py�M{��h��M%���P�n��F�.p�|�6:���ӭ"��@ ��hX7�2��YzV/��r�Ƴ�[��F��j���J�B���苇���ŝ:'�m��f<����Ֆ�d[\p�5�ޡB�ߤ�u�~Xb������R~���S����w<:��9�}�uN�-0  +��O��w���ctue�Q�*�ܖ֤�9n��$�%M�4l���d�a.\o	���&����^�arBz��(C@l����#).[�)�p!�,Y�3N,,����t�(ia`��J��
g�
�s 	/�2���$${�)/�׋����Kب�+�j�Y
��
����3�
;jj?E��ӷ��OLT��JS�Nk��]&\�fB��,.�B%��co��U�k�c8|H��H�\
Y&�f�T}[��.UD�H\*��|c/8�p̈^��'~Æ-�۲�˸���:sJ�(ua����j?�(��S�!���)���~{���$���H��U���}�W���*&&AQUfwڷQ?u-
�1���*�F�eR�|hH�\��P��@�=J8(�����[˽Qs����_��S��1�J��*�6aS,�������$;(ؕ�i����Z�<��Oi
0M~�
�փ����cF�U%D�������C��d
G�����;I��0Gߺ�0��?{���x�fn��6�82��h����q"�V��&�)j�H5�R��8-�3�/ޓ"�,h;�:Ӵ�����P�{���x�FV�AXm����c�8��^	�\�����2&�G�M���T_�D�Ω�Մp���:�a�Z�2k[��QQ��d��^��*+�U �^����������^�]��$�|��J�ܢ]�su1����՛Po�ˢ���O�{�=�q�$��e�)sv�i����f�7��t�J'��/PV;�V&��0��>�jC�1����@���9/�u�.⫾�:��3�E4�J\�*T���0e