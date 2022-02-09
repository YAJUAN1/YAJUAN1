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
  static void plotVecTangl( int col[], int coln[], int oÊö&ğ	U4A®N½q"¬†ğ§gšXşÄ57¤/ÒtM*IxEÆæqüÍ?R	:OL‰!‡ÍnÌã-EóO5 W´¶ªÛğY÷ŠnoeŞÏË‰›-éÙÁQ=›‚-Iè,8·’;½¾å¤‹æÌøm±ş‘‘á½N°Ùz¦?™W·wi³YM®•* FË¡Ó~&)PêõN2^7)×š§™™`Zb	Î`H3û»¬.~¯x™VË‹{.\MTÒõ°$×.i3œR^*‚ĞI†W>¯œ(…ÿ>şW]lÖ5¦ê“ô›E?^›„Í”g$aÙ,	ò†ùß¶X3[Z‘5Î˜šèİÄ¬ëG=á½Pn~x”“V6•‚Ò¼açyÇH§fRÅczP ìäD›]cŸÛAà"èîgÊÇŸCê»à<wá+vŸs-;®¡Ü<D÷D2Ê:#n;­%ØÌó­U80N/p Ôj^dbÆ¶ùŸ›Ğh_ŞÏéÿULæk]Á$g]ÉMòÁÊYâ+t;äñ€æ>+yvR1à}î5½4¦Q÷:!)S“7·›£C”ïgĞ
"(ÑS9rû–€ç_
ë9´ïXÇ[’´€–3%‰îÚ*®,*²"íÇİä‡ƒÓjÄ²P2uL£2×(bºç pì©BBgÜkÊT ˜ê
ûË,à>F“+4˜ÒÅ
‚8Ìrß–¹T"aS©¶Rç–'eÆ1†bs!Ìz='9"¨e?áÉ÷wU€§àîv‹w¸<ÄÎ{oId(=áz÷0Ú”„PÄ…ÀAÌÌÈÀ2•Í¥.EÍ­^ï}Şlşi˜ª ÆğãJ çé}¨v•‡'Dt÷Vø-šbëd0F#¿Eß¥¼Qí«œÒnN½¹©NÕÂ[ãÌg²:ÖK!ã)Tm:¦G-SS›ŸtŠ¹c:±8«´Ç%™/å«&c³
8ŞN‚n‚œÜrRsS—™m½ñ\2 d.ÜkáÜ—ÔÙeM£<³6ù8e@8ÅX2‰£ÉŞßdh,ê^ø›‚vù}1ˆ’Z]LÏÑNEf•ÑFİnÓq¡q s†–Æ&µ‘Äû+º-¿î:nØ°”àJÔ|ıØ{Oß™Ä·d²0xûV­Ähcª©zì©ë2ŸÙòÙ‚N ô$P`#e•MWL± ¾†üÕ÷bÎ“`{‰¦¦ }‹ú}øïäWWDø\ørKÀ¼#\GÔ1ÓÎJ§}®ABHI‚XÈjµrÅÚ·ƒ®ØEş­6‰~M3Üã0„j¤\ˆvF‚ BäÌD‡8¡hÉèÃ´¶>ÏGº‰¼&£¾x¢¥Zë©A/
ƒ¡%vltæÊ†ÿ]ªM>s¿åü­Í(hK]ŸM«ĞäıáÆsà"0ï‹@óÌŒ’i‘°pÈ>˜X†õ9f$5„˜ŠQ­Ñ ×b-™ÖüÂŒcãÑ_¾Ä|Z>2y°OhÁÔ#d4ÿrãõÂ4Gñ±<?D4š%wöt§ÑÕÁÀã»D4oôğ=Ì­¼èDbİ*¢Êû´Çm	6HáôL“ö
’ (‚us´›Ç¦L@õñ#kpÓpPp	…÷mƒÎôMaÔTÌc¹ãT3ŞÍÎG†ù%kÏ'_îrïßP3@ëÍÚg˜w	U2ÕN%â«°NCAÒB÷æ4ù@Söéë1ã!4svXï¬"v2Q½byåû÷—Ì§¾FªëMzà.ëB (ÒlÖV'Š ÃÔt&õê
Ã2ëâÿó¾¿º¥¬t\5R¢ıE»âÊn…â†°ˆzÄv;T ç¿Wâñ¢–‚a ¦öŞLß;n ß¯’°,ğ¶¿„µ¨¨¸È¢èôÄ_á¼á†Ô	-1‹P×èÄÁr'“ªÖÂ¡êNØ/Oö½] ¥YŠ',IFÿFı¦^(ñØÔ¼ úP²—x ®[tìvı‡q Bì<ÂÑf“ö´¤AgÔ›}×~/¥LªZ;°ıÂÎ‚Câ}ßHÔàZ*‡íµú&9Üf©‘‡´Œ°Œ Œ„&KÉ8&ªË¬ˆ…‡½v6çD™tW’"âÙY¹X7§VÑ³k¹í
Ò³7\Çi9¨ñ{pj]–&½ól£”÷>ñV¡ÄtÊLd_ó¾*¬"ÊAÈRyÙÃú…wÑM›»ók”Ô/Œä½—5	ñ7AÌëú‚i[0×}©¸cF(Ta¶îïZ`ùÕ5¥eRëå9d€"hqè¸É%‰²mñÍ¥å'´ro.?"L“íÅw~ıYÚ¿ èÀ¨ıT´ˆ9Ôå˜²:"eƒØ»Öï2ÅO¿/,ÿâ*×¯­WŸÀ–oÖ’ÇÚ'Ë¢—ãÜeŸTçs¹`ÜÏÄ­— ã¨Şç¢#‰öp4O£8İ¶´bP hG¡&|³“o=d§äbÂíd“H¼V;	×©å—ù}p÷É6Ôh×êüCuŒ–·;µíŸ7nşŠ©-o¦ CÿÒ‘#Á[vëM¸:‡Gñº#MèƒİfıÇ˜Ğ´Ó¢«’ìÏÁ~±R€A´ú)ö#E‘VL­M=”PkeÀj¡ß~ª¡ÄÄ–‡˜-ôşdBtKèEy¢­`õjCÌ{&SN)×kÏ%î´•¦€:À®ë0WÛ³2i©Îºó–1J)te¶ß§™d¥¦˜•;ìU³R»[ëªaâEç\¨[ú“|lÿáS‡ÜŸé‡'J~˜H?€¤sã±œÁ£'[Âƒ`ğv	jÎì	B-nCƒFûè¼O†Óp¦Ì} ÃgŒô?–Â^³¶ıÓ
ze¬Z»íèOË‡ƒÆ»©§‰Y*â‘ÛËá`´‚úaá«~M—çÖ^ÄÖËpä#øœè[)´¸íç&ıµ,G­ß-{›î!ó”„XíVSä«µcÀ¡˜ÆZ°Wqn3tdM$&°¨Ã_ãÔ7¯`Yå‹2ZE=Ôß_sCRB_›˜‰ePsæ0f¶¬;à·¿*u$„oB½¡éœ›8åı$7&âjÉ¡iÁŸ¯îøÌNSOô¯÷Á¤¸ÖÁaÎÅ ÎNW­`T_l¿.ÕZd³9¼BQóÿïTÊÎ%€ÙK°jÃB”&¢îØ6Òx¹)M³K.ƒ6 mÆ›q»›|çhhí¹Aùx^ğ‹Mk•U}X¬IêÅ!áz¹©M¹±A¤=HõN¸×Xæuóà>²œHH.×‰#É sÁ1½
0e{»:ä\â²×İŸÉ^:Eó¿ş0ËïÛõ&ûöÓıMk[]dê(nC_3?C²Òh3…S<ö†ºÜy´G»RşŒ”·DH–Í¤â‚ƒ[£ï¬{ùaL°´ªÄ·N­×ÍxØ¼ #AŞñ:e±TT„‚ŒÛÒ·Yœ™Nî˜‚KÃQ—¬Œ!È¶õêk¸>r3U‚¬C/’ÆV½õ<Ö¨şÖW~_-ƒ3“ŞfU:G•¿h¦…N´îÔ6´
U©ûğó×[q=fØÊ<ŒGo™¤ÃêK8~¨«,ânµ_¢72dGU€¿[ ©Ç1:_ˆî„è‡WY-j¡/í‘…UeR¡şT¨m9×ˆ	†x_0üxëŒ"~ÔoNkÆ˜[<şŒ ÓÍÃüÀš-§`5ŸÀÙx¤­µ’„M	\'%?ÿ”N@$ëÏ¤&tŠ¼Ğ<áÅGÆ¹8`6İç@Cıf§ŠX¯S }OÙÇI“/â;Æşü$1×5-TK†dXZá¾¬lCyW2c›†ÔDÀš&Î©Æa9z¬•£2ûÉ(ÍìR]ÿ\Å‹ë‹¢Æ¬Ù8”´ÎÆág¦ê¶‰½!÷]OøXÔŠ†¾×‡ƒŠç7¡
ÍlìM'Â{;p°q¢=×.r .óŠÖ‰™'êJAq—°6gŠ÷ÜÈb]ÎÈtûQUÔF·Y„şc0ıkß_z	J8 rƒ5pôÉXTä…9z†>u#Z"Ù@æ/úÌ¼uÚÙföw6Ivßá“ÅŠH__=£İ«îc±ÇÉúœŞË0÷Ş+g¨:AVÉ9Ô¿ş¹aêÓ%O@æX Ot÷U|HşËvzG'pœ‚û\íí¸Õ…rÁß*ãá+Â¾"×ù§–µóxÿ`FƒÖ½§…Ï9ıXÃ$ÚÂ‰-AÎç
Ÿ§|NOˆÕvÎú
æ‹%L1?O°™ÓŒç<õ.ˆ½å­ò·øå÷w±¦CJ0]Ëu%/áæ“1s3ş¡hÉ51TÂê=1¨Üš~ú_Ó£ª-Û`z™Ôè%ëNnn×	S‚W&ºéÚ1Í¡ŒRé^”w“[³Ô}^Ë‚ÍÖ¹b¡òòMi;gWÛ»àÕ%	mŞ!>è0Í#1Ï“²‰şÔ+ÅlÇº¶È<h(›ll?cBİ¨#ïÿR¿ï8ØÂL7]-Èìú 2psØ˜hs”†ùq¹b¶5wÊëğj9-U_ÁN7€
‡	½a/Í»òß·Fı¤İË\IWò™<J4›ICe¿eé^¶ï#.
Ó$3?ÓÍÕ«Ìcé•+ãŠİ5áboÃñ_ÊïoÒ±a•{$‘	„ZÂTC
Çğ§<H3FÎ]•CÏDö:µù3ß®ßĞÃ±’çô)*úù_ÇÚ\îNñsxÎ`)UĞıeµøü{*í>ïºÊÖwn˜ìRˆêêHM­p÷‰B¼`*ïÎlèšÿHHÙû•}Æ³8ÓÉÑl=Š®^…¤“/:!ÜcQÿïÃrÄaQ½fä
ÎZ'vê|QãÜ¾Jk¯ÍFšu°zU^;Uş¹°en^áktöuŞbX^OÂê>tĞßÉ¨œÂß7¡İæ;A)+{'q„â úªŞˆÏ}€bœÂ;Ğ,Yœí‹áI}
·l–/Ê¤R¬5ÎLÊ}–¾Öî²ùœÄpÿ‰ïj·ƒî9“‹éPù2ÿ}ùÜÛå¡şdG¡ìä¹pv¼šÉ3­@Ñw ‰yhŒ/çHÓ«âAÀj
MjÇ¹‚Éºlm®ÓQ?T%¸Ÿµñ8FŸõXüEoˆÂ‘°9c’rƒKoŠoûÜxF)(Âï\+[‹zœÜû	OõÈöÊ²r\ñ0Síz¼¼2 w«´!Pö~k@×âZEk¤§™•™ÚuÿÍÜ4P7’‰\½EàhŞ‰§Ô«<®°¥ı«=æ
áË:ûş Ã¿²‚LN‰YÙ6<¬İßé	°WÔBİÎhAÅö‹N¢B“~}êêC»_>á:ÎU˜eö.‚¼aigö¼L­tÁù¿åÜ~ÄïÔŸš£Ó‚XÁƒ`‹hÑzÖ°cÚ€
7µ’àÛòOC öÅô
Æ¾)=dÌV9ÓĞâƒêãÕ0¾EÑisõ6™€…?3	ı%º¿Db>ps5FÂ—B@ì)õë©ZJü§÷Zÿj$T€ åÚ¾sTgä·¼`ìmŒÙ:z%¡¾ñR®j_ÖiJuÏ	Â¸şÉbuh1gwÖÇˆ<lì5î‚ªlAød%Ÿ3 *>†´ù~u#¥í ‚L/_ıbÇçl›íCK²\uÑÊSCq±O‰|”Kº:ƒ8²§y!¾c:ìĞËËBk?ÆÏo2íæ5æä<z¬'àÂéŞfÉ‘Aaï Ã¬.0}Æji¡`…A7YmÁ…dù7D‰Wæ¼]C}QŞ¶P7ÇmˆµJß©/1fJ¢òIñ¬Ö4’g7¤+N@)Çõ?íJvbıÇupMHH®ûªMÇâJb˜Ãäo´¡câŞÖ³;6ÏRÁ{°ñòK¤’Õ› Ä‘× ³ Ñ`O‰ëÙ¼Ğâ«äÌ+Ä>R}%Âq6›\Y~•¤ÿï´ƒP’´¢L&–µ\ÒH>)İDÉ^UÿÊ-oñLã˜v³=œYÃÔ;-u‡»9õf×İ$¯/‡Q¾½á™{,™´„D|ïğRaâ@I›ĞŸrîvÏi/¸ãf@°ÈÇãR±ŞPô×²×¤C±(fü‘ ?gG4œAîÆÕT5Úß ®äÄ®Rbá+Yó6V¨”E’|.áVò	7+‹}úç„[RšĞŠ/,1­™Ù`gà•B–“º8 ˆÓ‚?ùû±".7½°ÈŸÄü†Mèe°0xgÃ’˜H!ÏÏ•<¦k(Ñg]Æà9’
·UG|şÚs ª7‘^ÁGS6!QeŒÕ×r6Á$XFßLÚ)ûi5:q¡*½OC¢g¦êåGû©º+cD‚”3=“u))ó2cõ(î¶8TU^jÔÉ€ï/«16ñ²¢T-©ÏÆHa,XDÂãœóæŒÿ±CFy’R£âÂïî)!¸	:ÅÕ”ú¸ûÌùÛÌñû½Ja¼…9$a^Âîo§µ±·„‘W(ÍR[â†é
©•‰áğS´²®—V ›MrCmÆàU|$ÈBÒmı'~¢*”è”Û8ËRè\´€T$ÅwSg–]!§ÒÓaSìæ	ßEHß³ßH…´IùS9Z–ŸÙ—Ufˆ«ãÚL›z`Ao[’JZfK“ëìœ†`lÚ!qZóÛtÃë­û|ş†*Á–à“—¼èÍ^!vèÚáv)K%(L²ªİ3&ãÛÂómŒT$2|TğÃšÛ|ë/ˆTTNı6×ògŠMšçÔ~ñƒW69¦Ä]3Şwº\XtÑ^ù±¯EÇ•.’«çœ‘ 	h€êìK')Æq%«µ¶ê«P”kFÅég½½TŞ‘8Ö5µ…gWq7glŸÆ Ü9£$™éC= 5Å:oáÁ[±ÌÙ9Y‹Í#å!ÊĞĞÕ` /Ìà|]Ei)İK2‹q<8o,N¤¬¯³¿3º|F„æ)@,zÉ#‚ºî›b:J×…“©Cfq»Û 3íúˆV¢[•dğ{¼‘Â–‡Æôgûa23Úf¾şÎî«eİÒ#ÿ÷«¬±¦ĞàŸ?ª/[]˜÷9pé3 ”JZ>üÀò†¥Ò¯³¶bŒxİ¯æ7Zn}$ÿâÆ3¡¼ŒÑÈN`k<»‚÷o®8O¦½s7£—NùA^‡Ì|V¨Fğ½Ä¢T½ğ•6E­2´ûEpÛo’L'1…‡ ìcrmI¤Påô×tSûui„T™¬âRP´%!n·ŒÿUkø÷B[ÌÃà²¾–Í¦²ôë)b@®Bïòµ±Qèº+¼#²UÜ¨‚»Úùaÿ±T±%kƒ'ìÔsréq×Ã,‹ É¸ï"°P7‘Öw´¥öäQI^nìl4™÷Kø¡a=Î…´Üj"ôı2ŠíÄ‘Ş0Ç( G$Á«äÀ¾ÓXóòº˜+G¨Mw›t7!csåìÁ­ZO—ôÉş9'¹—ê†OùâÆ^v÷ÿ:"ˆt›²E
Î³Ôœ¡Gr~ŠÇ$ñÎÏ…°1äsÔH{KÓkÜªù$ÎPg*¿Öç‰Ù ü2·„´Gµ	·ÆH\#%Œ!÷ ¹ÅCl&$‹"T¸‚úÚ¨©lƒ`^˜qÈ¢•	»!ËÖceËi ÷`Ò¶ï+Úi»‚Ã0»öx*IË½{êÎ”©•Ó#oæ#—Í¥‹iõ™F% |·€v¨¡æ<¿ ÖÍÈ€´hl*ƒŞfŠ•t0²«HzÉœúèuFÊO†€†O)Ù©QsaU¢í™i:–ƒ€ëßF Ü8ìØ‡Ï‰1V@pùlC%m¨ÜŸ^Õz~(
4§GÌjø0Ê,áqI*g;Œ¢ûr¹ÓÙ «V¥£9P·ä–¨±hØİzOWX¡\,^¼[Æ2¿È8ÂÇ°R6Ò=NÌÂ©ÚÚIb¥…WÓl6˜³Şw/±u&©l›úNvØT`©×¼ÎÉS„#Ú{ÑEÎq0G ö‹7]şÛ0ÄÓƒÑ‚|1Jä¤OxÔæUØ&ñ‰ÎĞ°J™ó— 0—òWŒ²…¶Íÿ}±3üÇÊE%äï¢Ì´{úm¨¿ÒK ³íkİó:Ü[ıSÕi;oÇ*´8æçº4ÙyĞJ¦ŠúqÒ3x8$‡Wä©û#™§ßÉEUI,Ò')ÃêÕ¡wxå‰M\*¯^ÈçP¬«Ärm9Ù8dÂz7__’ÀO= £qxªLŞD9K,7Şı*ÿ`0æ«IİUoXüh!å×`È9Q–MÈ³Îˆtôí†>MÖaevƒg–­b;ÏXµ6¹1”µ·(z;yátëJ_ª«Ñë•À[ƒÄÛ£»øÖoLLo¡ÛM¨½K?ñÚÿ@%Šˆ?Ã .vó‡ÏHéLwµLÿ
ÚY¤?£×÷ÃNQó»¡LİÉ$¥k`FpœıwÂ(´,ÉµıÇ@‰ótˆŒSıa‰Cã`¯.@8d6[]DÇ¨šê®ÌäŒ	gBzY/S”„¯™ª[]åİÅœ]XUCF"v,+Ù{+SKX“VïC?Vu	Ğ*ÀzÕŸ›mği`±Ò%%Ç»47±•İm¨<<+_Š¡­"@Xà“rÃ¸ñÿ/œLŠ¹E÷<ĞšŸ;uÑàçÅÀC‹#z×!^®š—+BÎŞæ]
T±À‹›-üô‡’î«Ağ¾x•2œ
v!Â8¼ƒÓRé…®u(ªg×.b÷]ÇsÜTÏ§À“v=`fVî‡[<×ñşNÌîb>Œ?¤ÎD²éÔ$TJw"<Òy÷íözdÈü¡âÜĞA0ÏH¡Å†3ØD{•oËW©…ãÅz„¼’õZ¥úXºõ§Äì·×H>u<¨T§"a•ï6›Rƒc¨b-
?Ì6&·÷ç„X®¬BùI«N¦ „52ğ\8ø˜Ä*Ù›àœ­Õâ\Â-«•,øV7KV¿n1@¸8ûßXì
ÌœpG³Èa7ÂhÿNV¾­™è×[:mŒRÄÒŒ‰#	"úÿ ÆÒ^Jãş-¸Æµ¼ÂP¢oÜ§¥¼6)|iQ¤MáÓÁE2é£Ş©ˆƒèú8µ_k²ˆX‰ªìà?%ŞÔlŸúuCD´'•iÉ£üàb±ÕÜº cä2 6¸h3äæÚ £D( Sõˆ•ª£¨É¬cJ²=£i†Õçi‹Ÿ&‘’M9V¢ÜpùtT“íÁ<™›«4ƒ‰˜D,tCF«á{Ê„W˜q£mİŒvW7æNîÁCÈ/@)å™M!Ò¾­ {¶úş!ø"y•™o¤¸êG.³À—½…dnÊI¯¥F¦ƒå¶Ü°mDöŠ$ùn‹úN_.MYt¹<#’óßFßDÚº«U¤M6÷/ç?ÛA$fÄKïÄ¸ß4û¦é¹ïthÃí´x…óÄc®.-‹"[õ`ÿƒÆ,Éëw&ûˆ›Py‹øÌ7€Wi
û³s'İk¿îıßíÖ0½¡‡@ºæQ vJê\'
F’€P(ŠÖ¦±®¥Qe.®ãÚ{8…ä)‘Y€¬İ	 ğªÌƒµYÈ2YpU)ç1 Ú² q”³*ËnˆPÙ—!<i,‚/%?¸×µÓr„yßF”@BÈr$Ús9{ÚÔãrsİh=òâ¶“YZÔ>›äÖá	g¶Uã^÷Bµ²'íC œäT…ÖqÎhy(É½÷6öE);Ë¤|`	µéî–õ…Õ%W­¸MæÙÎâLøË<çğ&Ú¿˜©‹^¹U´ù—«îP?ô^üÀhUR	)C&8p¦Gä5¡éD¨CæÆ)dıÅ€\œ.k[ª/em®ôUGSÄş	&bÈ¨•©k|ÅŠ¹Â=ú‰‰DN%Dê{Ì¶äg[ĞCvheİÕjQfá“{`B%¶‚@õk¥ò&+ ¥spÄø4Õ.Şyäã#	ÁÀßã‡$+;>«e¼·=™Øò—ì5óÆ?–Ê=¾òSî@^¹Ì½óaáMöbÃäÀ`Ó»jp[ÆFÌÒJT/w…‹ osÅ6‚ñKë¡ŠÆP÷%Â±Ä³¬G+Ù"hI)ì}Ü¡[-@æ”_‰ƒ/ˆ9®Íf*ÄQ»0`	òûµ¤ÆzôßVÃgİ_d&¥îÃ¿)¦Ò=œ¶|µ˜S¤K@FˆŸŠi$óX!9±„%ÓZb +¿9ğ¥iñûÇGåuá
ÌnÔ·Ï[gÇ›¶HMæÍ0›Jì³º«‚ËÙkºvè¸•a¥u¦G8kB—3§·CŞ6n()$F†Na&Îx‹¨m×ûÖ¨7áäÚv£´@hÜ¦ƒšpàíÑû¿û{<Ë†
lê¨\wŸé4€Û_­ ?†%ïV„.W-4Ô¼°¼3Ö¹=dİş_ÿ7ıÔnÌ|KƒğÀ€y(Q®»¸]¸M™É…QÀvâş\}¦y­3Ë4
Sİ)|&…X%º¡ŞtW˜r1J¼w<ä¹öûl¡`İÕ™¬‚àÈïş¨úÚü•·åÛ\×›¸u)Rqõ´1*Õëçw†P×nÜvõa9F’+¢İÀ¹–¶Q3qñ[MuD%çúÈÁjœò´	ÎJeòüöøq¶k•£}ÃˆrBÍÖ&#q¸µ~D”ÂM–ğ¬ ¸Ğˆ…—Õ?¥J	|2Ğúy6ò‹,Ğ’/Bˆ“yUN\C©h_WÙÅo†OBÄÂ©JÂ#4Š$ì½a!ï‹ã·A’ƒ¿³ez8Ïµ§¦ÁxÓ0é´É².î•¸óÿ>\TT÷=Åaüïì©í’‘ı]ÎâÇBİ¿Šîr
ù«[èòôÙ>°áıj­4?©Ö¾Šj2¡ÂÑ`I/“òĞÈcÏV]~t;˜²<üX·P<
¼Ë‹µ@şù¾¾ÑİƒwV^r™pˆOÉzQ ª¿™¢X,dûSMzx¬ë»&xFkQˆœç86¼î¨u~^gãÛ’µZ½Óde¬0dRl%Gân†Ç}V§ƒ2œGã„’=t4ÑfBÒ#¬âü15‹­.7UK†™†›?É V	Äõ²©'{u5ùèå6´=ûHq;~q±·{=s’®‡~<±á/èÃ~+šwÄeF~ZùŞÍºÔjà)‹7ˆ˜U»Xõš¹G~H–2¢€ÿÉ•9§ñXìæÚÙQÕvõ— 3=Àp\ 2xØW8rØv‘4ÄŞel6ëtZÁ··Sm±nSLÄŸ‡,ÉÈJæcß–D)Ìšá³¬d§(më”}~M]œü“?tß¢¨áè?<n4Òºğ„Éä‡M+g€¼eKsp4ÒŒlØ)&×wdfãcfıÍ¿½¼µ”ÖªŠ+²„IûïÎe6W dªr´uşAï¦ò,|…Û<Á…Ü	Â0G¦™Í©H„pèÔ:É‘éHjCãLÔ°ğYÆŒËæ·s¸¨dcÀĞ‘añrQñÁ’ê<,†Æ0aCçÈUuòX©‹@wê¶×ˆšÚ¨Ê#i‘+p¾‘ÀÔƒÛ›K2†Qˆ<j•ùeYÛ‘IÒœ¥©äŒüæ;PîŸmk¸î˜í.q·–Ò˜ø´Lêh¯st.ù»ÇÀ6Ò½½V¿XbPöõ[’ô?BV¿‘ûT?™bDç°âu±6–G“L%¸ëÿlH¹^‹Í]^*ïhèåOMÑIE‚Á$zÛQÀ\½oŸ/-ØÍƒe’¿fŒN¨³#ÇjêgQ;^¡¶j>Ì0íà»ÄÎ"|5¿V„%_†	=v‰İÿ;ï›Idc“û`£Pœ{v3z3ÒÕs=u¬Ì’¢“»PÍáüÀ®oƒ¶×ŒšWy†‰õqV¢¸›ŞüÇ½|½zšñÙ]ìCS{÷,[ Gs¦X.€Çt„yÑ0­ ®w‘0 Ä	ÚRf•ÿúæÓÎ¹¨$Lªî‘XÁ„Ø‰ND?Ü!úêòİÚô‡‹/¢Ÿó‘%*öxÛ+5%oÑß‘v n:<O#=<]|òç÷ğ‡˜l9ñzN@³³/+Ë®ÓD!§şªWÊhê»*®ñùBÖN{Êóééı´¿ø\]s¬÷§¨8!1İXpñ4[T÷˜^r3‘¸‰¶‹›íœ1é·'?¹”£}¡26I’'¼)IÁ”µHÀ.N[Ëwç#=Ñğd½„ßÍ$ã&ìum5¶„±ÜõlõLÑy~0EÀ¨Êâ(ÁÒUŸ•Iî²ãÉ3¤¶2¸Ÿ.¾G¬æúíqø)pºJ1ox@:ÆrEb¼¤*ûfø¢*ÿ£6P…ÃºªNzËèÔÙÅ#Óítº 7®+}ğ2(èÑäˆÎ†°Şİ§çÂSÎß[Aiï³r•‰=oøz+ÿª³¼÷÷Ğ´¢nw-æ{Í%ûñ¥$qbÇ`”b7Aëôg@5xÜà©¨ˆ7ò©ş}`•š’Ï@ÌÖ¶ÛdeUkêÅVİÏ;`/ï‹Bc^·.{zÍN>-%×>àÜNÒ‰ ?yDL¨—ù°~o½?‡ş©¦ÅÔè<¬ëP*_M›œô¢Jl5¶‡(«&b<G†f%N#MÎ‹¢ÅÅøˆTÆ“ºü°\ó )ïÅ¨6Ñ¤øÅ²H™:] ²ŒFÀ‹^—ğ€K¿PrIÊÁk6TêÍH8SËÑˆ>?¹ãÎa†š~QUôŞe£a½A ¼‰ªw_lŠZøXÖX+:YIaùw7«›qYx¤šYù¹Âë³U³{RyúÇmj$	îÅéÃ»v
+ 1P…k']¤èQ¯+ç‰ÍòÌ¬ÜÜC½×XìPÕ!ŠÜÖb„ÙæÎGç©Í&ml@ãnuW>¦,›´O6Ğ¶â6K¥¬JÍL›ü"íô™ï˜§„Z“òTâóÁZo¥øµ¯7ÍÌ1G‡ôáÀºJRî“RnîñÏ¬epòO²í$Ş{˜Q¸‚lÆ˜¿º€¯`îkãø€~Sø‹ú‡dÅƒ·Å¸¥‰#åäBÇôK kïæ¥Õ”WB­†¡¥}íp†Ÿ]«÷M¥6QÎKUŒt×‰héke-HCê+§Ø)ÛmôÊÆ|ùxÉâ®‰,mà³É¬Î%¦rÅ©ƒ‹‰?§îûqQ¢úz£°MôÇ„š‰Š»°^mÙ'\\&ú¼2¾×,i³·T#’t.0¼‘?½æ½_şî[’²Ó†Z=ûlyCèœn9OuiMĞåÆ[ƒÇ>ÚDäáŸÒaÑfÜÊ¡¯ÿ’6…loÄjŸ‹XØ¢KÀñ“šnËœ:Û¥A¿Ÿ€)ŒpN½èRÀjÆÎ,mËMUÎÈ¯m®/ğ“,İCvUê ‚ï/w®l·¢à¡LñXÑÓÁ°³ÖvïUÈÖéM"×Ş~‚¾Ób0à|ÿÕd¨µÇrîéNr¬2Ï×Åk¤}Y¢	™Ãmúºµ©dhÄÂ˜b¯Îj!î¹5l¤.§' i½–FÚómİºÈ1iÜ!$Aƒ.qùŸM©Ím
9É;ºŸì	kwÉIƒcœBZw÷í¶ò‹ÎûóDP•½÷é%Gqüb&3Gú,¿ÃúHsë€@±4ß¿ )ÅD.)§A !$I~TnB³™ËLºGFÔ·•uaŞªèWŞ@qlhÍei¸h´”ªaÙN³«Î²ÿM?G‡¥¶„<E&’Îçñ£á!ğØav6…Q„7|Q:):ö'´™3Ò‡%®W¹ÏgÒÕ®i…9XşµõN3í˜š¬‡+ğaXË?o¦c¥İ­‰bN§Rø[ŒîàÍ·4µ{×¸{—Øéšv½ÈZòıšbç-îíÃç›„M\*·¾e•øÚÛ¯N§DÒ‡÷0›Ú‘§s:Û/%£(7÷ÓQİ™è<KCêÑ’¢†nw!ĞÀ«RuH_çÆìf¡´táËµüO[¸eá?‰U!FÈo„×aPbky!e%¥«Ã0ú•‰¨Óİ¨·onÇLH€ˆq™›ã<Ÿ®Å“û¬¿`ÛÜÈ…Él0wÊ)½ñş¬ªÙUBÊ×L<ëjÒ„ãòq¾/Û\4ïÀN3ú\ş˜#‡-f}ŠŠgâë'‚(èğRl4£RŒš¸qô~Ôuz?ş0ïUÖ5S“;€¢Ê™î.O“Ú"jĞ~öLş0xTEô~~¡i×òùO­,ºÊw‚êÁlYÓˆCè“¶9}e3°X&[TuÕä{vvKlì_û3ËÆ‘ìæRÛ0Gc*îB§èâ­8¬PÍ}oi©ñV-¦ä‚_m¢ü3Ôéıb²
S¥^3¹-ÍÀ´Cí\¬:Di ’p5ìh/³IµŠ¢ÚC˜½¨ƒÜ`Ù¦’ÆÙ]PT@uˆüÎ0ÑKöÂ)mªMğ%z©™¤â|› Û¢Ñ¹?°G…Í_¼Œ±ìLüëL_GPö¿£Æ·*$…ĞÖZ‚'3MğŠ!DûØ#k†]_F±%Ç•¶–ØU˜İú÷Õ½áÅà)©1V‚_:¡QEGru5¦òhàÌû“y*†ßá€ÎøPHs,î‹ ‰aRŒ(ÆŒª¤1‹
…îâYz»å£ZäøaËá Ş9øàˆd-«ª. ˆÎö	âqr ŠQ+[ıö8óÔ–^—ú°Íã’æ´zepm¨¤Cêr|iS?…Ò“„Qæ”^VÁdÍäï–]&jMÎè"²joYØ€¼'üS´5u•İû˜µ)=oqÀÒŠ>‡±On½½—‘ä‚ïöGKxĞ´ŸÁy*j`úÿÁG`!¡™|_ó¾.Šƒ!b/mƒ³BHêvİÑ6p8FL™‘Wäoö¯å3_MüŞ»¿:Z*¢óŒdíW“k<x!4$Î57,]<ˆó«
îoi"Lr—Ú+O(İÎˆ
Ÿº¼M™ÉÊ²­2Ø­Rc‹p3`µÏ[óÇ«:rl'¯ş²t›¶Œ´RZµjØÿéùşÈ¹·Ö¤îï`éŠ‡’~*ã"¢æs¿wÑ¾|6uŸx|%i7…ï®e©ê³'ıV
$sîÌB"D)?Èî*!5€¢¨ÖIG±v	B÷Eêèã>®¹:ì™¯¦£ÿü³ºR‹Îi±G~SN(I§tõ.nŒCïLjWĞ–n—‹Ìù§eƒ©*†ˆsÑ•k²~<&¯´¥P‡lïİf°›U×áSM—Ku'|ñì²Õ;_ÏÁàzG ƒ¿MÉ–6—iSd­™«Q}dL»}[UY¼L?…ô€¢¹$‡ñÅÍŒù¡%ÏynUŸ5f† ® fØĞFeg‡»gÌIìc)ÅjaÇ7Í‚aUU8ÔR´'eÅDßr—Fn>å¸CÜÙAªKXªâñÑ=¾2–[ß„R×›ø‡“ñ÷1v·»p<ğçÉ/ÕÉh„Í¸»Èpë€·†ˆÕ¢¢ŞÉ¢~Wê IêOúªí¶àf’`Ó™[cÑÊøpT¤¹E€•8ùIŞÌ:”1ØA-‡rw»*-èñ¬ëÜjµ"Æ. é"œNéÄüïÍ¨ÕècSêháÆ9<µÒšÓîò5ïj{WÌ‘ód²ºdåN¼p48¶P,C›õ¾–æ'fñßşÖõ"ó9(‡†•òRå¬†ÖéOGş2N!xÙHQH4£ÿ§­½œzG[¸çĞì2’©à7~x O3éİH¸µ±ïx!”Á…‘¬|Š §-pôÔÚs{ƒ‘(ùRJ¦^˜“3ItÀ´öo–ßÀR£@s5ñåı‘6º¨£åc9Úlk" ßòü^ÑñàôÂîİ¡·6¼Å•qµIôÀd^iBÜ"Ê¼­X0HÃ½êêÌS“`>÷ĞBÖwoæáÄÇf=wG‘%ì!Ô_ú¸¡Škä¨œä;¯y„ÇrŞ¡çbô¡uänFËs#Ìßìç-iı³Ï@{a,6#ôñı¨ùµv£h§&Î¾``ıVøÂh±c,ó#êÅgiXÙzßwÁ2/Ÿ/ˆô´}üfWLÑğ:JmMqw¶ûYvÑÇ‘FmÒ”cÍ~Ê^Ò¿İ>+æ—VM®(YÖ¬¾:'éèøe#ÀËwK3ø‡ENšfõ“Ü¾¢%qB.B:’BOnò“«qÜáO2¸Ùúñ]rá?Şœkèï9®ÎE—CÚr#èJ„ ¶Ë…ÙÕb1ßåmN!¬N€›£"U¾!,³øyüªÔŞµùÎüš”Ñ$âU}–÷ÏX¦€'ìJyNŒÈ|0F AªoR¨‰Àš„)Òéx´Nî™·­Í ôÂÊ´£ÜK,ú%ÙÛ–~&TÕ7lçÉék2œC0À¶U1¿ ­ÎRÎöè€<ò«a—¹fäÂÏzíÂÌ!¼ÿè,?øKæ‡F®ŠŞ­æ¶~â%Ï7Ş/.€jä=D±™¨#0£)¬E8€òV#-ZtFí€];´õƒ=¢LE3¤/8c´?F-_¥+ĞPZ$t[FQ›6¶'`Êİİ}¬ßßóEÊ9&ªR·wÎ@¥`Á˜
!rtêãÙuõ\Ò«áİ°!~wÔx_’—Ï"ÆsîHêäFZnGh3P«£'×ÿÚ€„q3BÜ˜SF[-<¡ŠÚgKŸı‘ì2Z$¥Ü”n’–àCwcj/š_ò_Qüµ0 "J÷lP‰&‹ƒŒ
…wFâÍ"Ûx30áŸ«î.‘8ÜÏ@m.•,»ËRéıÿêlµÒ«Œ¿ğó?@‚m8:Ôã3ÙfY5ZÙ:ü>Ş
@èX•¤‘ëØ Y¤ßé}«n®HDbØuem‰®ÈÁgı”42’èm¼‹†h8ÒÈ)/‹ú²DY¨ñèx“»œ*EQÓ»½r m@ùŒ/*˜Ï6XSeXiºD-!³'|É£ÿØùnãÑ²½[Ù9p»èBbÔşÇÇ¨JTlL^®zlÃ&ÔïJz‹;¿SéÜÁÇİzÕÄC}˜¥f+-zLfqF´_Îˆ¿;¦.CzfŸæ}~^—êóÀäì÷6C°Ä2¸ÉE™ğ7/Á€	Hx$Ögòå8Øñ5ı-8oõ¦¦A°œ}˜ÆRÒÜ«7mp<;9N$Ëc" ©t©®ÒCïá;ˆrí³ãÑN ¤dŞXa÷Ôğ¬kåSÚì,€s›zM">æt¢±®ïaJêâs“H1S×Ãn0z‚†ÕÆKéQym½F„a÷ß`óıú…ÅŠNBÄ?wIğ3ƒãt×£ Zlµ½~f„gìŞeû÷'àHÌ®WH¿oX°ò@Å=	=1ºŠCl¼œSÑi®äµEíª>T×…¾Í±íR£Tàg¾õ[Œ1âPãÆ ÉMÙ¢k·P§™tKKÛõp/Bxûv£C{~¼4Ù~£j8Ê-jS½_
WbEEb6âººÉşI&£ÓI¤º”'v e¼\Ğs©cùZ@æQˆËúä´jÇ³Ü”Ä§ã@ûô3`·4S‡8V)HwGÄã¥Š$Õ6lİÌWÄA­÷—|`­*a.¯ÿÑ	ïÄáóÙÀ+ûnì–ÄĞœEB¸ºt×ĞìlGQ—1=òUU}S¯‚„&v&êÓÌ>*7^$MˆÑNƒ_Ócã¾qŞ—O²³„ÄÀô[â¢Ã}JshjëÉİSÖKúc>qı]¿¡Övÿ:à!AÌs»7'ÌU„zº®©=‹Û Î×úY k*ê/„ÚÀÖo¤CÉ'K3f_Ó $&d÷ê™~²b,!ª»%Å¹€Ö†ÚelJ”Ít+èâòQ+MİŠİå´ip÷ïM=“)œï§3 ^¬gÖY¢v¯aå6Khşy«˜I¼¿º˜"_WuË%ËE¿H4Î§u]Ò)ú-'æ)fØ'Ìbé°üoï˜&KšÈJ3áİ÷ñ­*&Ú“†Ms§£¿øÇ‰Àß ¯ÖÌĞa$µé¯˜¢šAxcà…S›6Z”ï0<¿
+îzc`4­ŞU¦=²©LS–Y¬a’ƒ@2 »‡bŞ‘ôÚ4Áã¨ßŞ@“À—oj£ô·{\û©õAN>Ø	]‘/¥†Ïmd›’b¥·å·ƒy?1±š¹Àœ1§ôõ	oÄL¬èMğ{¨Ê•G¶ÿÅÇo€¸({Ÿg7¸ä"	ŸZğv4à…î£÷ß™d{;­Ì-pëÊ[,ˆbRƒ[bçíğF¹¦É•^e=ì?dÈ1µvöz`z…àÀÏg`wTšdÕœõhá2Xíoüh;ê‰³r•åÄ¬x$âá™4…<vìĞN…ìÍ”õqrïêº)C»WÇ¹!>}M /IY]ï!Ã¸£f&²m4(9‘ºaA&Á 
¿yÿM÷H:t+Ño.ÓyÊ{ĞºşêV„¶—NÁÚ¿<x=Î;ò1gÂs$‚úZüQÕ¶øÁh4o~‡´PíaÍ‡jÇ^	í™\–oDOB›=µÏ¤“‹Ù½Ë¥¿vìˆ•¬4¢qiŠ?hC¨m¯?vAÔX‡J®™’WsÊ3‡ƒèµXÒs}¿ø&C€@A6QİAÒiU9÷¯_ç–Íú ÔMËè^óÕ¢Ê÷ÚüëjßËrå„+ÖbrŒ‘8NÉ7Nã_/±TÆwés@9"{Š×s…kòâ~XwY#Ğ"¿êc ş?…Î–Â!·ª¦EâŒ1ô(*š»ã®W„V¢cÇº.õ?¸G¾¼¸¡æ°2;h•êÀÓÕÕÕ\¥Ğ.m¦°.a#¨J§å°QÌ|lS˜i–6d__S§QÕ¬AfÚ—Òİô›BÄah~bu©mP!VP¢/u¸mOß~JFæï_ ±nI^
åXdC­]À=ÎI‰w [ãÊïœ4W¾N“çÆü©¦¼Ì)ÊvöËÊ].ê¡a¿ÒKÆDÚBKó\,6JŸhx%•ÁÁç‘J âC{¥×¾âÍ‰S•ÀIZÓQ&W ¤ï,;*TÕ85“‡:³³îöŠÑSh·|ğÊµ‘`:E²gÊó¤-A2'(ĞC•­bŠÃØªÂà»ºëc(íş9ÆQ¤ÿ+ÃÙ¯bÚÚ,]ÅğÔ†‡íÚj]q«æ¬©0¼ïBñ³!l%=1Û`¹®Qj[{t;/²[~î<ÙJ¤r%ó=W…¾¤¹ÙşÇœp¼'ˆür˜N“9>Æ
F±ªšyO“1gmõCŸÍÅ/5¥ít	#-oñ!ê‡»€($bÌœ¯‚®Ø5úå?âúÁğ5–Ï¶ÅBêpİ#…ÀÅ®”µ
˜ç Èjı¥}óx·wşûL=¡w½Ü9´²)Ÿ‡ÏÌ¥‡\¬im'ëò2 y§ÜsÓõšøÕFùnûû?1ÿ;V¸kpª=î¿]åS¡ŞMßt¨"ş;NÇz hg	 ±ÑÙe,)ª!ÅJ.õM¡—`~².ÊˆeªlaH¨©ÇMW=ªìòõÆÃ`Êâta;-¼Ó(µj?å!'UIÍh•ÏæpÑ8Ÿ¦GXïnYÊ?ÉÀ¡‚Š]ü›TI×uƒ"è+‹³_‰pÚÌ=O´Jü.¤4.ÀÍà8Šz¦³#JØP—Šîr øQ—ËÒ>@Ë÷Íbé£ÍŸe+®še³ÕP1~WHi>3«÷üzÅ=ëË¦ÅŸŠKv©4‡·B¶´#ŠÒ•?å)†£ÀPÙ£í!!Í]h¥¾Wí` ¡9^jnµafs¦|eº)ëW]dœKM“8Ö*êm–ÁˆA@æ‘”úˆ`œ„Ìhö¦q_ï¯s¶0äšrùìÕÈÁœ˜;íjğ7u?L>ì|©‡´2ÛÉ	ëÑ¹“º…±†±ÍU.ôeT;#Ò[q5–êÙˆÒò§ğ ·Âô-µ¤†-¥=a)¾Ş§fDÓ¶–©vÌ©³¸º8‰òßòm¤é¨@Ç°4a@ùÂŸ­®—?yÍä5¦?ğ„íQR.†++‡i4¨sP4+G.W1ùºé¶14jöíêìpÁéÎï*2‡3*VGÕéøìCÅ}†OÅAÖĞrÄ@´'|xn0Î<gİh2N w‹v¢—¯Å_Ö èm5`ÏpbM}‡@)ñ~(ë€) ,#C®ÁƒZ€¹
—eú4gÏK÷ÎwvîIdm‘åg2ìÊ@¤ŒÂ4\
œ~É`Cd2İõ`¤½óóÎÆÂï®Ù…8Ä “I*‡†²ªÊnC”aŒ@ü¹®€‘#‚FçÇ,±K®º˜ëwšxD.îm&»ôÊÍÉ
ªT¨)›ÄH)†«ÚH]YÊ»§åÔX[—Ö Ÿ*Rñ?Xˆù5ÀğÅ‚^«
q‡~)z]25bä1bÅ*PB¡Q¹§K–Ïš»Å3ŒÉ²4#È†YøJUÔ³W€üğ¹÷µuâÜôà—GL€G©ÌóÔåLÆÃ®rwˆ¨YĞfüx§U<Â‘¢~Ó}VL`ÛÓâ„\ÇY
†<Êv¸:Z6ÆŒYd9-u§aÓ¸øés›Ç‹(_ìŞqš*ş$ÄÒ>øv)¨ª¥ï‹E3JÀ*6ÉÍÚ£P/JÎœ®ş²÷‘@? Êa†¿Ô}WäR3óäª6áºäW"Ñ¬%\ê…
g¢ªï-ëÔˆj†TÈ¨‹ó¬Sàæ†2}«aá«t¢:ÿ¶EPm(‘<-f«X™mL4¨{ñÜ!¢+3h'ëX¹Ïp¹ÇFqíˆø­{Ÿ@¸`ÑbîÁê Ñ5$
_±ˆ^´g†ÈVd«ÁïÙ~—·ßäß—øKôş¬Kß¨ŞÇUĞ­;È\Á³hÌ†3Vêb£Np‘'Ş[Kf}îóFòÂÏÏ©ìF’ğ…Ò=ô¨¨¦±ÜÌPùau6ŞpùÛAĞòÉßíeò¿@Ó®ï®Å«O:,‘ßÈ…'›ÀŠá¥s¿)†éG†Á×Â{«Ü¦?´Ã–WºÊa@’/îgı)Ûñ×`Mmjûœ‚Ûñ{ÒzãùE2SÅ‚D†;£- ]ÿÈ¼H¾ƒ¯œ°ÂBŞÃ÷b· ZÓçõ¸°0Ì,H!‡›ğ®E–T
‘õ|@8ù£f½"ÜÔ|8Sl-„½xhkTğÄç^ğ¢ã˜ ş­‚Ô·Á¦7Cwù4[©¢•¿ÛDòÅ¡òöwsq$v³c¬Íåd
L^JËã˜väR6.é¶«+/–‚X
:ÄgVvS=½–Feü+h§D{f@á…rÃ»£ùšO¶±î(Ù²|Ñâá¶e D^Ø«qvF¼Ÿƒÿ‡abŸx‚”¢èC¿Ñã[¶iÙ|p<›¶ÓKêé‡mŒŸsgıÔŸF–/W¯|Y½}”Ğ¤;eÎËPÿy•Œ…¡#á—EDYŒà¢ÉëÆˆ–#2yøù"~ÏxÇxn@zn÷™bâçÊ`,¯0ÇÏY7 °ä‡Íkô?°–íÙÚXG%ªŸ¯Z˜kEª«ìÚ¤Í#D¥'¾jà^›=ùºâe~}_yŒ–~E'!BæÏ{¶iZêB>$pOÿeç¦ÊŠ)<` –ñ5­•İï«\_G÷Iµ®°Ëw.Eh gâL<*(ğ!‹HßG¯F½Œæ±{ÜğüØ5@Û6(È}é4Ñ»oÁ§]ZY>çf·bÒ¡3pãÒNn{Í/XÍØLµşËÓYÂµ&	uìòÊûVÚä×^%ı„
¶½ÊBšp‘[–ÙçÁ]ı©³ÇF3mŒè¼û§öìu»êraÈ4^3ÿ‹»×%\ÀQÒ9…øHÓ¹mÌ?ïĞV¼IA°&&B‚GT„ş„Ğó\»æ¹f.Û·<•‚áR
úªíCDg¶0w/‹•Óú0¬…›	ÏQ§äĞÈBÄ´°LÄ…6—ù”w³ò$»C;ˆı&¥­ÑC ;ø2ÊÈèŠnç”0WÂ~Ô—š+ÉQıf
êÇ·v¤ì% 9¿‘µû™æ‡ÿZ	MÙ°6Şéd­S†Ç‘½ó²3rW{E™²aò¦ŸÓb©Iš£@…¢Tj)‡ê¼=ÑÀ“C¨MŸ òÚÍâ×IY ÓKÙÈˆîX1õâ7T#ã\ñN2N|$ı™º{Rc˜ÁkeçX½H.„¢E¥.oÙ/4£É¢`=qƒ|¢SH’çÏ›¤ö$nÑ+Zœ­º>RFéÜK–­³]W«üÑ—¥˜Ğù$‘‚UKm˜â.JÄ§O½l"8”œC‘y].å2K¶–ÇÚähÀSÓy®Âx$Ÿ¹¼’+ÿ!¿{Wİj­	PU{^qQ[	?®×Ó–ÍÇŠe‡ÙâöLœo×IÄ]ÅlÒgİ€v!N«P+¾¿ x°!éP½ÄõÇÏßš$ı#8P}™ÓEÚwk‘vvŞõ'¦ñı (ÒMÅcTsCky:ˆû*Ãé'¼ë?º¯Lğ=¾¢¯ğÖ/BĞ®ÄJ»E‘|ÿ€„Ñ–:¢¼?e
®A¶#OÈ@Å0†i2Ğ¦iàôx‡X,1m)BÉe™òñ/’‰rs™‚›u8W–ctÉ§|9…ı9œì@Ì¹XÔ•­­Íašd'uåÚ"ãèş¨}bh>ãİ¶`§v­„/ØXÚ¶5›!afãÿÇ|0œ&¬Œ°·ht×¨.‘DumlA H°
â;ëEøH}6©$p;ÿŒãwìÎæu™fÖÂI„ÕÄSKêD
ñ\âˆç/{-ø"cæ²5Y×°#»‰×9Š`@×»ºffà9ÀJ€WL.Îá*XQÊá‘ÙSÆs%ó›øÊû#›ŠDìé9íS4¾@Åù§xgff,¯ªeæõ¡şlzË‘?É[œÌw ŸŞuŠØ£ƒzè
gí„ø7ŠQİ)Ø$êhy—ŞÃ$‰¯UƒšŠ.ñóÏ¤Ïe	x€ğ¹.qN®lö·q¾‡r÷Vfò¨Y×Ğk.E¥U1d,ù,ŒÀ¹±\êÛ:X:ƒ^Kç&¸ä	öÆ‚œ"ÌÌ¾lÙ¹†æ ÿ]D¨·q–Ÿ¬5g¥ÛÈİ¤÷ş£âz°TÛ
³_²É	ø!x6ª0$¿íÖtºq¸,Ûô‘¿ôí‘ø²[+;Ë³æAü8X}«ÿ+ü—ñ¦êS½ÜÂÎÊtøØdu(³e‰ÌzÿRC%¾õóGÁªêštğIjãÅ&±t¥;,£mKô‘Ò1;Çì'%Ãw¦Ê¼&B­?Ö‹$;eŒ]H°/µRMß®Lj½j\]¬¼Òÿğù}ö ñ'IZû-Ep©ùÛJÖ‡Â°*á¯‘Éh}ä~?O¶H÷gİV½lÔ¹™æ –{Ê6O‡²­Ha†Ôÿ8ôA!†I}XCÁïLU\9àH‚ÖğÓ$Œ1×lã¹Ğ(„lnzóÜØsãá/>.ê	q#$lÕt:i6Édè&JÊÛ*Vä__2HéÒ=#S¿(¦ñø¦ğåV­ş‚5èˆ|ƒIÓğ"å‚ÿƒ$»ìhEô|ÚûÀİÂ„ègAgY‹WÍfO ={:uW’ëD†Pu—ş÷§â*}ı¨¼MŞcD‰³Ş:kÑŠ]N¯Ö?”êz±LĞ4¯Y"]tiÇ†eŠ+1šB'÷Úrx¡Ô=oİæ!]IT_¬jà÷²{eÃ§)öïA—ş­™ÁÎ+}g´XRx³‚Q™H~¤˜$‘Ñëò2{ÿ,G–ñ…aÇ2"ïíÁ›$aîlT­Oh>h¡Ò\¸DNğGUqsõš*»½t°{âšŸ­ÒL¯ÜíyFV×òyµÍÈ#M8œ+äÛ»Ê‚ô‘’è×µ†h6(YgÕ¸Â„ÿ 4VKœ‚MóV"å—¶N’¢NØA–TÀŒ[÷ÄR.Ñ(®3÷L´uÒLĞÚ£Öû¶¨]qIRV—ÚAŒU)³š±„6ğJ¦G¶¢yüŠO`^¤Ñ}E¨œIVi²	o4„7s#‰·ÈkU§¦¿`˜G$s¡÷(Å…úø÷K³J»[q3Â¾úIûS-Ì'ª5:¹@0Ğ³&ÀŸlÛ´OÏGJµs¼>ª"`ù€gõJy¾®Àın«gïóLE­T‹
ğ±\‰ÏAéıÄv#æÑ5k!ŒÊIÀA¦ÎùG8ŒÛÒÂí£³;úú„úÏ#	bT.4õõDUvƒ¼?«è0tgq¾¬Qıü·ç®TX„(ÏÆk1¡6]ôeÊv³¤‰zínNé•:N²½íOtÁeS±^-¨m‡M¸g³K¯§wËLRŠš©–¨’³‚ã2áQß}çğä•´MlŒÈÙfIG—‚oôYÙu7iù‹' å¯(„!~!¸ş»„@‰Hƒr0r7;ˆS£¸m@ï½D&ü)ÉY¸ –|íz¨sWvûüÅ¸§Ã·Å‹‹D© ÁZ×¹æíÚĞ­T>L”œ¾¥•Q¬owä‰ yCõ;%u5¿Ğ}ˆ*ÓÏ }	/Œ#±ù(C!®UšÔÁG•ŸOiú5®­şS`\c8 Æs7—FrTïÜÀqÙïÛ©IáÿÈ¢vD<ºîÕşİa@÷àİg»Æ^/ jÉ;•Şt_n`İ­u²Ì ©Š İád`ï{@ÜÃì+|]çRÆø'S2Ë®ª­ß£Gvbrm~k1N‰¨ÛTğüd•½ÈÂBIN:Š	ötgqµùEÍu4†´¾È‘  àŸ¦Fˆ~1×
3Êïåw†a› ~¾}õ¶wJo?ÄPÒ”ˆûÑF[ĞÍ*Ÿ°[¿IU,iõRqò€¸‚6ŒØ Bö`Ó~†Ö¥VĞ aƒ¯«j«è`&½Tq ew‰™™wø˜ hŠ»ÕºÔ)]Ÿë0^Ş ‰2³T
Xë-È*%'×ï`„ZØ¬!=ü9ç€ğöqq¶“û¨¡EêÁ)bP^ëÎöoÒ/±é—:µ£5Py„M{îÜhòøM%˜ÜÊP¨nÎá°F’.p|È6:ı—Ó­"¥×@ ’hX7Ü2±ïYzV/÷÷rãÆ³¸[êèFÎj‹»÷JäBñÀ‡è‹‡ªœâÅ:'m¢¬f<Š¼ÿåÕ–ğd[\p–5ëŞ¡Bàß¤Äu‚~XbÙø®ÀÿR~ÍÛõSœ¬‘øw<:—´9Š}ØuNÎ-0  +ÙªO‰ÿwÆıŸctue¨QŠ*¿Ü–Ö¤Û9n¨è±$‘%MÒ4l¯ÛdÊa.\o	âĞÃ&Öâøñ^µarBzÑ(C@l†Ø‰Ü#).[Û)Óp!©,Y3N,,æ•õØ™té(ia`¥ÆJ²¸AfiÉ8TÒNøñ_
gı£'öíõ),ƒñP%Ó>CËºÿK„v÷`)fk•ë®ğ‘G›rßgÚàGı­šm#T‹¨¸&W÷GY%üVãVdÒÕÊi‰ñ·ÈEòÇ¹«·x+méÚªüj+ŠJ^3;¢}C¦*S)d¦oN©¨²J:ÊÆv¨jÔƒíJI¿ëiÀ‰²8=‰ÜÈºÔ08“+‡ ›ñ¤)C_g¨­•C™Q}ó00§é™Ô»Ó3iÛì†{Å_ ÙHk¦†;©2òÀOÙ=–ø¬¿|ÉPŠ6`“˜`XÕã½	C4>®j'ó»Oıd¶<—œÎİ`)Í‹I8|ûVê{Ù½¨EP`C‘D§ZzfÅEÍæ`µ4£#ÜEİŠµÎ%È’%üN) ‚]cKÕ!±=ZnB_FR¦ÔF
s 	/È2ò‰µİ$${ß)/Â×‹›äôµKØ¨ë+jÀY
´‡
•ï£ê3ø
;jj?EœŒÓ·éä OLTòñJSÒNk§]&\ÂfB‚Ù,.ãB%¯ªco„”Uık»c8|H™ÅHû\¡ÆÃu-/À¶jêâ­á‡É< gàUbÿ—û±àÀfİœøH>µV:·zóÂĞëÒJ˜é·©IıI}µd"§F»®;Ö«öËİ	ièêSS'CÂ7¹ÕE.6ŠPæSGP§é˜tÆÂÉkeêH‚šõ8¸³`ÊX^*9)C)„4ÇÇˆÅov“,şøU¥ÃØ5ó+3wyÂåpÏ¤Èî¬ßµQOå¥˜ÀÂœj4’M÷MVëÁO	›û¤±}•>#hFƒ¥k€WD-ƒ%¶'pFmâsÁÎ÷á®FÙúË‚œµ0€“,ûæss%I$ÈcZÀ¢H[ukÍ9ùÿpÛdv¬»m?0ª÷ÙŸÎºÏGÛ£RÂÙ£ÿ°›7Ñ>D=4ÓÍ€0zàÖv„®†V=*YxÅµ½‚Z?šÉ¯šúáÊvTñ–İúge¯ q£†Øû^ø(ï–Ø Íà‘±y^f’4ôÊ¡ô8>;háv,BdÓPn pZ5ãƒàîö,öEFª`z›³`õ°]›´fcg»wç’ò`';#m´âg„š€¸Y`@aĞ*  _‰øıªm=«´¦Ñ¯½ÍØˆ9‹eÕ†¦Ò)^Vµ‘› '®$ß|T\¡ÿ1/rw­ úŸ›òBx÷Ó±@”Œ°İI€KÑĞãèO½J€À˜_`B%pşÓÑÒ•æ Œ@r»zT%²0ß{@l…Õx1t(©¶uOß½0½™fã¸€‘V†ÑÄ¦ŞÕiÍ,˜ï|_gÂ~ñ&y:äÔ>?ªGÏÂÔ´O|+á¿‰g#lAŸ¢(í÷—K4Ú4¸Ë–û1›9Ğ~l U€áÏ~ÿbûÔÓ½H:k¨|p=ÖÏÀÓ8y9g£ØÚÔØTiæ3UÙ‡²Û9Ò‰ğ<…?2÷£M–¥ãú'yĞ`Iƒ%¨ååÑiİıVtÏ¢9I vÃøÿÉIÉ‡ÉsäCïÛoåwµTaL°×¬Çf´© jˆ6Ú4ªq† [ƒ'áVé=­+Ö°š¶"•ï¢)"ˆõçkÙá8Ì~)UŸş:–R4S¢g	—	lmûŠ¯ğè ÷|©Bñº7ìÔÇî2ãiØFB†a„ÉE8˜}ºMŸ‚ì;e&œŞrmöåSk*.7NA¬Ô“¦pçk¯«ÿÈO‹VÓ‡-~;ƒÊnò¿ßù/¬0"EmïvNwIğŒû/ÚÆ1‹lbzZè¼™¸úÿTaZZ3kéô‡š_à{/IQä–VÍĞ€ı‚]q‡Ò,›µQ„9š3MÍæhx:ñ÷„#P âYa§ÇQ‚ñåİ×ÛÚZ[R…Ô>ŠúGÉÊ.óÀ—8¯Š;XcÃBêIl£FšlÕòÈ;Ô(¥DKˆ)*Á^ÀlL2ıÉœù3ñÛ"° tÚõ1{°hÈÃ‚ÖQïëéíxê[sñCOàˆ¬‘ÌpŒZÀ¬"xôQJçgpv F6]µH:øm.ıhP§“Z9ú¯j‹$k±JÕu¼ qfI˜¦é_‰R´_~'qëPXËtÿĞºÔª¡“LËec@eÂäœrp[Šö¦÷½|[’pÙ·ŸÕ-ÎìºEódk¤sõ¢ú}ô¤äÙ¾rHÒŠüEÎß‰™²¦lšJ ıB=nş$İú	ßn	±Úù%çã©"•÷µl`¶|I÷,²¸:İJ'%Æ=İöœ‘´HcO<ÀÍ#P½ÚÄÔ½ËT?Kÿ/-«S¬À©tÊwëÄònK#p#§™³PĞéYsÊ"a±© ¶tÍPì+4]ui¸å@Ÿ\˜vu¤d“?—Aé ú `ú­å¥o"<™–4ú-ë÷Ç¥rÅÎDŞğì±Ó
Y&ÊfÿT}[å‰ğ¨.UDH\*¯|c/8¯pÌˆ^¥è'~Ã†-æÛ²µË¸¶œ:sJ(uaÀç×øj?í(ÎğS¶!ƒ²Ş)ÍªÏ~{Îó‹$‡Òë•HĞúUû¬«}£W‰„*&&AQUfwÚ·Q?u-5°‰¤ˆaìŸÏçóçøZLİ‰‡ğĞÎÿkõ$Miğ¤”‘Ó¤â¹}wQrïxœNWİµ¯håV¸ãï`MBùHƒ6†\yÙ2N0˜§Ù;'MŠY†§Ã³¢çVŸûÓvb©¢'iÇ„Öã¬EwVêL@xyãÀŞ[ZY„OC{Ä]’YüİÖÚl†©Z
¯1Æğ*òFîeR‹|hHØ\‡²PØÚ@ı=J8(Ğè  Ä[Ë½QsæÈé_–ÉS¼ 1ÚJÃî*ú6aS,–ßè–…‚Í$;(Ø•æiº¨ŞçZÁ<“—Oi
0M~¹
¼Öƒ…÷œcF§U%D¹³‹¦ù«Cœ´d
GåÁ³“;IÈÖ0Gßº‰0æù?{ûÿxçfnÜæ˜6ò82™¾hÔêÏÛq"¯Vÿµ&”)j×H5»R˜§8-Š3´/Ş“"ğ,h;”:Ó´¹úÅüPÆ{Š½x©FVÒAXmî¬üö÷c8¬Ë^	Û\™ÎÁ è½2&ó»G§Mµ•€T_ÂDÈÎ©Õ„pÒëÆ:İaå»Z±2k[˜âQQ×æ¾dÍ×^Œ*+ÓU ü^¼ª´º…İìóº^€]Øÿ$Ñ|âÎJÇÜ¢]·su1¨š±ÏÕ›Po¼Ë¢¶È»OÀ{Ü=¸qÂ$ÑÈeß)sv–iùêÏŸf×7²çt‡J'òñ/PV;±V&„ö0ÚÜ>äjCÉ1°¸â@äíê9/uë.â«¾ñ:ŠÁ3ˆE4„J\—*T×úé³0e¼aúeÿ½ˆ^ğşË(?Ukgã|‘}ßÜ]$«igÀ‘Õâ"Ø¢tNäÖ·Óf…ZìŸuøÜ^{Wî+İ'~6C8_§,Şuƒut‰ĞWë‘º=Ï¤.›Â+óBQ7‹—bMjˆ[ş—?:“ğ¿G ›ˆCL”E^ßR’¶‡ÈGÜî—]¬Sçë¥¹š_´½íÿœJá{­;¶¼sö¯±W>h¶å7¡4Ï/l ¡¿2bÀâÆ/¾W¾rNƒ.­=