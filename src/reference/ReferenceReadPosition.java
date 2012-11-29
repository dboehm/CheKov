package reference;

import algorithm.CheKov;
import net.sf.picard.reference.ReferenceSequence;

public class ReferenceReadPosition {
	private byte[] nucleotidesInReference;
	private char refAllel;
	private int indexOfRefAllele;

	public ReferenceReadPosition(byte[] nucleotidesInReference,
			byte[] affectedRefAllel, int indexOfRefAllele) {
		this.setNucleotideInReference(nucleotidesInReference);
		if (affectedRefAllel.length == 1) {
			this.refAllel = (char) affectedRefAllel[0];
		}
		this.indexOfRefAllele = indexOfRefAllele;
	}

	public static ReferenceReadPosition getReferenceReadPositionInstance(
			String string, int posInContig, int flankingRevBp, int flankingForBp) {

		ReferenceSequence rf = CheKov.getIndexedFastaSequenceFile_Ref()
				.getSubsequenceAt(string, posInContig - flankingRevBp,
						posInContig + flankingForBp);
		// the constructor initializes the flanking sequences, and the
		// refAllel
		return new ReferenceReadPosition(rf.getBases(), CheKov
				.getIndexedFastaSequenceFile_Ref()
				.getSubsequenceAt(string, posInContig, posInContig).getBases(),
				flankingRevBp);

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
			s = s + String.format("%c", b);
		}
		return "[" + s + "]";
	}

	public char getRefAllel() {
		return refAllel;
	}

	public void setRefAllel(char refAllel) {
		this.refAllel = refAllel;
	}

	public int getIndexOfRefAllele() {
		return indexOfRefAllele;
	}

	public void setIndexOfRefAllele(int indexOfRefAllele) {
		this.indexOfRefAllele = indexOfRefAllele;
	}
}
