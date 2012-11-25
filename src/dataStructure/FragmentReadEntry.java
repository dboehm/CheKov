package dataStructure;

import java.awt.font.NumericShaper;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import reference.ReferenceReadPosition;

import algorithm.CheKov;
import algorithm.IntervalAbs;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

public class FragmentReadEntry extends ReadEntry {
	private static long readUnmappedCount = 0;
	private static int fragmentReadCount = 0;
	ArrayList<String> cigarTokens = new ArrayList<>();
	// this static variables are for manuell calculating the median
	static int[] effLengthCounterArray = new int[2000];
	static int[] rawLengthCounterArray = new int[2000];

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

	public ArrayList<String> getCigarTokens() {
		return cigarTokens;
	}

	public void setCigarTokens(ArrayList<String> cigarTokens) {
		this.cigarTokens = cigarTokens;
	}

	public static int[] getEffLengthCounterArray() {
		return effLengthCounterArray;
	}

	public static void setEffLengthCounterArray(int[] effLengthCounterArray) {
		FragmentReadEntry.effLengthCounterArray = effLengthCounterArray;
	}

	public static int[] getRawLengthCounterArray() {
		return rawLengthCounterArray;
	}

	public static void setRawLengthCounterArray(int[] rawLengthCounterArray) {
		FragmentReadEntry.rawLengthCounterArray = rawLengthCounterArray;
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
		} else {
			// get the offset for the absolute localization
			long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
					(short) (this.getSamRecord().getReferenceIndex() + 1))
					.getOffset();
			// System.out.println("OFFSET " + offset);
			long absAlStart = 0;
			long absAlEnd = 0;

			if (!this.getSamRecord().getReadNegativeStrandFlag()) {
				absAlStart = this.getSamRecord().getAlignmentStart() + offset;
				absAlEnd = this.getSamRecord().getAlignmentEnd() + offset;
			} else { // auch bei reversen Reads ist absAlStart < absAlEnd
				absAlStart = this.getSamRecord().getAlignmentStart() + offset;
				absAlEnd = this.getSamRecord().getAlignmentEnd() + offset;

			}
			// System.out.println("READ " + getSamRecord().getReadName() +
			// " IS NEG "
			// + this.getSamRecord().getReadNegativeStrandFlag() + " ALSTART "
			// + absAlStart + " ALEND " + absAlEnd);

			// sieht so aus, als müßte ich von dem Read ein IntervalAbs
			// erzeugen, dass ich floor() übergeben kann, um das richtige
			// IntervalAbs- Objekt zurückzubekommen.
			// I can not remember, why I have switched absAlEnd and absAlStart
			// here !!!
			IntervalAbs tempRead = new IntervalAbs((short) (getSamRecord()
					.getReferenceIndex() + 1), absAlEnd, absAlStart,
					(int) (absAlEnd - absAlStart), null);

			IntervalAbs floorInterval = CheKov.getIntervalTreeSet().floor(
					tempRead);
			if (floorInterval != null) {
				if (floorInterval.readHitsAbsInterval(tempRead.getEndAbs(),
						tempRead.getStartAbs())) {
					// System.out.println(floorInterval.getCoverage());
					return;
				}
			}
		}
	}

	@Override
	public void analyseQuality() {
		// CigarString is splitted using StringTokenizer
		String cigarString = getSamRecord().getCigarString();
		StringTokenizer st = new StringTokenizer(cigarString, "MIDNSHPX", true);

		while (st.hasMoreTokens())
			cigarTokens.add(st.nextToken());
		int posInRead = 0;
		int posInContig = 0;
		long posAbsInGenome = 0;
		int numberOfBasesAffected = 0;
		int index = 0;
		cigarTokens.get(index);
		System.out.print(this.getSamRecord().getReferenceName()+":"+this.getSamRecord().getAlignmentStart() + " : " + cigarTokens + ":");
		ReferenceReadPosition rrp = null;
		for (String s : cigarTokens) {
			switch (s) {
			case "H":
				this.setHardClippedBases(this.getHardClippedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				// the corresponding length of the tag is always one field
				// before the tag
				posInRead = posInRead
						+ Integer.parseInt(cigarTokens.get(index - 1));
				break;
			case "S":
				this.setSoftClippedBases(this.getSoftClippedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				// the corresponding length of the tag is always one field
				// before the tag
				posInRead = posInRead
						+ Integer.parseInt(cigarTokens.get(index - 1));
				break;

			case "D":
				this.setDeletedTaggedBases(this.getDeletedTaggedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				// the corresponding length of the tag is always one field
				// before the tag
				numberOfBasesAffected = Integer.parseInt(cigarTokens.get(index - 1));
				posInRead = posInRead
						+ numberOfBasesAffected;

				/*
				 * we want to check, if and what length a putative Homopolymer
				 * is 1. first identify the position in the read
				 */
				// we have to subtract -1 because a deleted position is not present
				posInContig = this.getSamRecord().getAlignmentStart()
						+ posInRead - numberOfBasesAffected;
				/*
				 * 2. identify the position in genome, Inteval AND/OR
				 * IntervalAbs 3. identify and check the missing nucleotide 4.
				 */
				rrp = ReferenceReadPosition
						.getReferenceReadPosition(this.getSamRecord()
								.getReferenceName(), posInContig);
				System.out.print(rrp.toString());

				break;
			case "I":
				this.setInsertedTaggedBases(this.getInsertedTaggedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				// the corresponding length of the tag is always one field
				// before the tag
				posInRead = posInRead
						+ Integer.parseInt(cigarTokens.get(index - 1));
				
				break;
			case "M":
				// the corresponding length of the tag is always one field
				// before the tag
				posInRead = posInRead
						+ Integer.parseInt(cigarTokens.get(index - 1));
				break;
			default:

			}
			index++;
//			System.out.print(posInRead + " ");
			// System.err.println(s);
		}
		System.out.println();
		// calculate the rawReadLenghth by adding softClipped and hardClippes
		// bases to effektive Readlength
		this.setRawReadLength(this.getSamRecord().getReadLength()
				+ this.getHardClippedBases() + this.getSoftClippedBases());
		// System.out.println(this.getSamRecord().getReadName()+ " "
		// +getRawReadLength());

		// store effektive Readlength in ReadEntry object
		this.setEffReadLength(this.getSamRecord().getReadLength());
		// sum all rawReadLength to calculate average in main()
		CheKov.setAvRawReadLength(CheKov.getAvRawReadLength()
				+ this.getRawReadLength());
		// sum all effReadLength to calculate average in main()
		CheKov.setAvEffReadLength(CheKov.getAvEffReadLength()
				+ this.getEffReadLength());
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

		// calculate median by count numbers in Arrays of length
		rawLengthCounterArray[this.getRawReadLength()]++;
		effLengthCounterArray[this.getEffReadLength()]++;

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
