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
	ArrayList<String> cigarTokens = new ArrayList<>();

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
			readUnmappedCount++;
			return;
		} else
			; // something useful to be done here

	}

	@Override
	public void analyseQuality() {
		// CigarString is splitted using StringTokenizer
		String cigarString = getSamRecord().getCigarString();
		StringTokenizer st = new StringTokenizer(cigarString, "MIDNSHPX", true);

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

			case "D":
				this.setDeletedTaggedBases(this.getDeletedTaggedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				break;
			case "I":
				this.setInsertedTaggedBases(this.getInsertedTaggedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				break;
			case "M":
				break;
			default:
				// System.err.println(s);
			}
			index++;
		}
		this.setRawReadLength(this.getSamRecord().getReadLength()
				+ this.getHardClippedBases() + this.getSoftClippedBases());

		// sum all Cigar-bases
		CheKov.setAllBases(CheKov.getAllBases() + this.getRawReadLength());
		// sum all Cigar-hardclipped bases
		CheKov.setAllHardClippedBases(CheKov.getAllHardClippedBases()
				+ this.getHardClippedBases());
		// sum all Cigar-softclipped bases
		CheKov.setAllSoftClippedBases(CheKov.getAllSoftClippedBases()
				+ this.getSoftClippedBases());
		// sum all Cigar-deleted tagged bases
		CheKov.setAllDeletedTaggedBases(CheKov.getAllDeletedTaggedBases()
				+ this.getDeletedTaggedBases());
		// sum all Cigar-inserted tagged bases
		CheKov.setAllInsertedTaggedBases(CheKov.getAllInsertedTaggedBases()
				+ this.getInsertedTaggedBases());
		// sum all Cigar-deleted tagged Reads
		if (this.getDeletedTaggedBases() != 0)
			CheKov.setAllDeletedTaggedReads(CheKov.getAllDeletedTaggedReads() + 1);
		// sum all Cigar-inserted tagged Reads
		if (this.getInsertedTaggedBases() != 0)
			CheKov.setAllInsertedTaggedReads(CheKov.getAllInsertedTaggedReads() + 1);
		// sum all either Cigar-deleted OR -inserted tagged Reads
		if (this.getDeletedTaggedBases() != 0
				|| this.getInsertedTaggedBases() != 0)
			CheKov.setAllEitherDeletedOrInsertedTaggedReads(CheKov
					.getAllEitherDeletedOrInsertedTaggedReads() + 1);

		// printUsefulDebuggingData();
	}

	public void printUsefulDebuggingData() {
		Cigar cigar = getSamRecord().getCigar();
		System.out
				.printf("%n%s%n",
						"=====================================================================");
		// ReadName
		System.out.printf("%-10s%-40s", "Readname: ", getSamRecord()
				.getReadName());
		// Koordinates
		System.out.printf("%-10s%-6s:%d-%d", "Interval: ", getSamRecord()
				.getReferenceName(), (getSamRecord().getAlignmentStart() - 1),
				(getSamRecord().getAlignmentEnd() + 1));
		// Sequence
		System.out.printf("%n%s", getSamRecord().getReadString());
		// QualityString
		System.out.printf("%n%s", getSamRecord().getBaseQualityString());
		// QualityByteArray
		byte[] quals = getSamRecord().getBaseQualities();
		System.out.printf("%n");
		for (byte qual : quals)
			System.out.printf("%s-", qual);
		System.out.printf("%n");
		// Alignmentblocks
		List<AlignmentBlock> alBlock = getSamRecord().getAlignmentBlocks();
		for (AlignmentBlock alBl : alBlock)
			System.out.printf("%d %d %d", alBl.getReferenceStart(),
					alBl.getReadStart(), alBl.getLength());
		// Cigar
		System.out.printf("%n%s %s %s %d", "Cigar-String: ", getSamRecord()
				.getCigarString(), "Cigar-Length: ", getSamRecord()
				.getCigarLength());
		List<CigarElement> cigars = cigar.getCigarElements();
		System.out.printf("%n%s", "Cigar-Operators: ");
		for (CigarElement ce : cigars) {
			System.out.printf("%s ", ce.getOperator());
		}
		System.out.printf("%n");
		System.out.printf("%s %d  %s %d  %s %d  %s %d  %s %d",
				"Cigar-ReadLength: ", cigar.getReadLength(),
				"Cigar-PaddedReferenceLength",
				cigar.getPaddedReferenceLength(), "cigar.numCigarElements: ",
				cigar.numCigarElements(), "cigar.getReferenceLength: ",
				cigar.getReferenceLength(), "Cigar.getReadLength(cigars): ",
				Cigar.getReadLength(cigars), "cigar.numCigarElements(): ",
				cigar.numCigarElements());
		// Calculate AllBases
		System.out.printf("%n%s %d", "getSamRecord().getReadLength: ", this
				.getSamRecord().getReadLength());
		// Calculate HardClippedBases
		String cigarString = getSamRecord().getCigarString();
		StringTokenizer st = new StringTokenizer(cigarString, "MIDNSHPX", true);
		ArrayList<String> cigarToken = new ArrayList<>();
		System.out.printf("%n");
		while (st.hasMoreTokens())
			cigarToken.add(st.nextToken());
		System.out.printf("%s %s", "Cigar-Token: ", cigarToken);
		System.out.printf("%n");

		System.out.printf("%s %d   %s %d   %s %d", "rawLength: ",
				this.getRawReadLength(), "HardClippedBases: ",
				this.getHardClippedBases(), "SoftClippedBases: ",
				this.getSoftClippedBases());
		System.out.printf("%n%s %d   %s %d   %s %d", "AllBases: ",
				CheKov.getAllBases(), "All HardClipped: ",
				CheKov.getAllHardClippedBases(), "AllSoftCipped: ",
				CheKov.getAllSoftClippedBases());
		System.out.printf("%n%s %d   %s %d", "Deleted bases: ",
				this.getDeletedTaggedBases(), "Inserted Bases: ",
				this.getInsertedTaggedBases());
	}
}
