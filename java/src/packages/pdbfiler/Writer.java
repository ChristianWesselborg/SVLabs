package pdbfiler;
import java.io.*;

public class Writer{

	public String model(int num){
		return "MODEL        " + Integer.toString(num);
	}
	

	public static void Writseq(String instring){
		
		File infile = new File(instring);

		String outdir = instring.replace(".pdb", "");
		System.out.println("outdir: " + outdir);
		
		try{

			try{
				Writer inc = new Writer();
				int modelnum = 1;
				String mod = inc.model(modelnum);
				System.out.println("mod: " + mod);

				if(mkdir(outdir)){
					System.out.println("Directory created: " + outdir);
				}

				BufferedReader br = new BufferedReader(new FileInputStreamReader(infile));
				String currentline = br.readLine;
				
				while(currentline != "-1"){
					if(currentline.contains(mod) == true){
						try{
							while(currentline.contains("ENDMDL") != true){
								currentline = br.readLine();
							
								System.out.println("Attempting file write");
								BufferedWriter bw = new BufferedWriter(new FileOutputStreamWriter(outdir + "\\" + mod + "\\.pdb"));
								
								bw.write(currentline);
								bw.newLine();
							}
						}
						finally{
							bw.close();
						}

						modelnum = modelnum++;
						mod = inc.model(modelnumm);
					}
				}
			}
			finally{
				br.close();
			}
		}
		catch(Exception e){
			System.out.println("Exception thrown: " + e);
			throw new RuntimeException("Failed Writing");
		}
	}
}