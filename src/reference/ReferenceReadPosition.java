package reference;

import java.io.File;
import java.io.FileNotFoundException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

public class ReferenceReadPosition {
	private byte[] nucleotidesInReference;

	public ReferenceReadPosition(byte[] nucleotidesInReference) {
		this.setNucleotideInReference(nucleotidesInReference);
	}

	public static ReferenceReadPosition getReferenceReadPosition(String string,
			int posInContig) {
		IndexedFastaSequenceFile ifsf = null;
		try {
			ifsf = new IndexedFastaSequenceFile(
					new File(
							"/home/dboehm/NewReference_2012_06_24/share/reference/genomes/NCBI_GRCh37/NCBI_GRCh37.fa"));
			if (!ifsf.isIndexed())
				System.exit(0);
			ReferenceSequence rf = ifsf.getSubsequenceAt(string, posInContig-2,
					posInContig+2);
			return new ReferenceReadPosition(rf.getBases());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		} finally {
			// no close() method, so let the Garbage Collector do the job ! 
			ifsf = null;
		}

	}

	public byte[] getNucleotideInReference() {
		return nucleotidesInReference;
	}

	public void setNucleotideInReference(byte[] nucleotideInReference) {
		this.nucleotidesInReference = nucleotideInReference;
	}

	@Override
	public String toString() {
		String s = "";
		for (byte b : nucleotidesInReference) {
			s = s + String.format("%c", (char) b);
		}
		return "[" + s + "]";
	}

}
