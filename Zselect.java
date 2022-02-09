package java1;

public class Zselect {

	static void setz(int akarusa, int cld[]) {
		// akarusa与えられた時　その数に応じてcld[]の値を１にする
		if(akarusa==1) { cld[0]=1; return; }





	if(akarusa>1 && akarusa<cld.length-1) {
		double r=(double)akarusa/cld.length;
		for(int k=0; k<cld.length;k++) {
			if(Math.random()<=r) 	cld[k]=1;
		}
			return;
		}



		if(akarusa==cld.length) {
			for(int k=0; k<cld.length;k++) {
				cld[k]=1;
			}
			return;
		}

	}
	public static void main(String[] args) {
		// TODO 自動生成されたメソッド・スタブ
		int L=10;
		int kumo[][][]=new int[L][L][L];

		for(int i=0; i<L; i++) {
			for(int k=0;k<L;k++) {
				int akarusa=(int)Math.random()*L;
				setz(akarusa, kumo[i][k]);
			}
		}

	}

}
