package dataStructure;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import algorithm.CheKov;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

public class FragmentReadEntry extends ReadEntry {
	private static long readUnmappedCount = 0;
	private static int fragmentReadCount = 0;

	public FragmentReadEntry(SAMRecord samRecord) {
		this.setSamRecord(samRecord);
	}

	// Getter and Setter
	public SAMRecord getSamRecord() {
		return super.getSamRecord();
	}

	public void setSamRecord(SAMRecord samRecord) {
		super.setSamRecord(samRecord);
	}

	public static int getFragmentReadCount() {
		return fragmentReadCount;
	}

	public static void setFragmentReadCount(int fragmentReadCount) {
		FragmentReadEntry.fragmentReadCount = fragmentReadCount;
	}

	public static long getReadUnmappedCount() {
		return readUnmappedCount;
	}

	public static void setReadUnmappedCount(long readUnmappedCount) {
		FragmentReadEntry.readUnmappedCount = readUnmappedCount;
	}

	@Override
	public String toString() {
		return String.format("%s", getSamRecord().getReadName());
	}

	@Override
	public void analyzeCoverage() {
		// in this analysis, skip and count the reads which are flagged as
		// unmapped
		if (this.getSamRecord().getReadUnmappedFlag()) {
			System.out.println("READ UNMAPPED " + readUnmappedCount++);
			return;
		} else
			; // something useful to be done here

	}

	@Override
	public void analyseQuality() {
		// printUsefulDebuggingData();

		// CigarString is splitted using StringTokenizer
		String cigarString = getSamRecord().getCigarString();
		StringTokenizer st = new StringTokenizer(cigarString, "MIDNSHPX", true);
		ArrayList<String> cigarTokens = new ArrayList<>();
		while (st.hasMoreTokens())
			cigarTokens.add(st.nextToken());

		int index = 0;
		cigarTokens.get(index);
		for (String s : cigarTokens) {
			switch (s) {
			case "H":
				this.setHardClippedBases(this.getHardClippedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				break;
			case "S":
				this.setSoftClippedBases(this.getSoftClippedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				break;
			default:
			}
			index++;
		}
		this.setRawReadLength(this.getSamRecord().getReadLength()
				+ this.getHardClippedBases() + this.getSoftClippedBases());

		CheKov.setAllBases(CheKov.getAllBases() + this.getRawReadLength());
		CheKov.setAllHardClippedBases(CheKov.getAllHardClippedBases()
				+ this.getHardClippedBases());
		CheKov.setAllSoftClippedBases(CheKov.getAllSoftClippedBases()
				+ this.getSoftClippedBases());
//		System.out.println(cigarString + " " + cigarTokens + " rawLength: "
//				+ this.getRawReadLength() + " H:" + this.getHardClippedBases()
//				+ " S:" + this.getSoftClippedBases());
//		System.out.println(CheKov.getAllBases() + " bp, SoftCipped: "
//				+ CheKov.getAllSoftClippedBases() + " bp, HardClipped: "
//				+ CheKov.getAllHardClippedBases());
	}

	public void printUsefulDebuggingData() {
		Cigar cigar = getSamRecord().getCigar();
		System.out
				.println("=====================================================================");
		// ReadName
		System.out.println(getSamRecord().getReadName());
		// Koordinates
		System.out.println(getSamRecord().getReferenceName() + ":"
				+ (getSamRecord().getAlignmentStart() - 1) + "-"
				+ (getSamRecord().getAlignmentEnd() + 1));
		// Sequence

		System.out.println(getSamRecord().getReadString());

		// QualityString
		System.out.println(getSamRecord().getBaseQualityString());
		// QualityByteArray
		byte[] quals = getSamRecord().getBaseQualities();
		for (byte qual : quals)
			System.out.print(qual + " ");
		System.out.println();
		// Alignmentblocks
		List<AlignmentBlock> alBlock = getSamRecord().getAlignmentBlocks();
		for (AlignmentBlock alBl : alBlock)
			System.out.println(alBl.getReferenceStart() + ":"
					+ alBl.getReadStart() + "-" + alBl.getLength());
		// Cigar
		System.out.println("Cigar-length: " + getSamRecord().getCigarLength()
				+ "\t" + getSamRecord().getCigarString());
		List<CigarElement> cigars = cigar.getCigarElements();
		for (CigarElement ce : cigars) {
			System.out.print(ce.getOperator() + "" + ce.getOperator());
		}
		System.out.println();
		System.out.println(cigar.getReadLength() + " : "
				+ cigar.getPaddedReferenceLength() + " : "
				+ cigar.numCigarElements() + " : " + cigar.getReferenceLength()
				+ " : " + Cigar.getReadLength(cigars) + " : "
				+ cigar.numCigarElements());
		// Calculate AllBases
		System.out
				.println("ReadLength: " + this.getSamRecord().getReadLength());
		// Calculate HardClippedBases
		String cigarString = getSamRecord().getCigarString();
		StringTokenizer st = new StringTokenizer(cigarString, "MIDNSHPX", true);
		ArrayList<String> cigarToken = new ArrayList<>();
		while (st.hasMoreTokens())
			cigarToken.add(st.nextToken());
		System.out.println(cigarToken);
	}

}
