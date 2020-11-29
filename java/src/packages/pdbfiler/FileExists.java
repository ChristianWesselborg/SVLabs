package pdbfiler;
import java.nio.file.*;
//Checks if a passed file argument (as a String) exists in the file structure

public class FileExists{
	public boolean exist(String in){
		Path p = Paths.get(in);

		if(Files.exists(p) == true){
			return true;
		}
		else{
			System.out.println("File does not exist: check input");
			return false;
		}
	}
}