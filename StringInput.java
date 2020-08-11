package pdbfiler;
import java.io.*;

public class StringInput{
	public String userInput(){
		String out = new String();
		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

			try{
				System.out.println("Enter file path ex: 'c:\\User\\Documents\\expdbfile.pdb'");
				out = br.readLine();
		
				return out;
			}
			finally{
				br.close();
			}
		}
		catch(Exception e){
			System.out.println("Exception thrown: " + e);
			throw new RuntimeException("Failed Input");
		}
	}
}
	