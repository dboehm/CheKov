package dataStructure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;
import java.util.TreeSet;

import datatypes.AlterationType;

import reference.ReferenceReadPosition;

import algorithm.CheKov;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

public class FragmentReadEntry extends ReadEntry {
	private static long readUnmappedCount = 0;
	private static int fragmentReadCount = 0;
	private static int onTargetReadCount = 0;
	private static int offTargetReadCount = 0;
	ArrayList<String> cigarTokens = new ArrayList<>();
	// this static variables are for manuell calculating the median
	static int[] effLengthCounterArray = new int[2000];
	static int[] rawLengthCounterArray = new int[2000];
	private static long startTime = 0;
	private static long endTime = 0;

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

	public static int getOnTargetReadCount() {
		return onTargetReadCount;
	}

	public static void setOnTargetReadCount(int onTargetReadCount) {
		FragmentReadEntry.onTargetReadCount = onTargetReadCount;
	}

	public static int getOffTargetReadCount() {
		return offTargetReadCount;
	}

	public static void setOffTargetReadCount(int offTargetReadCount) {
		FragmentReadEntry.offTargetReadCount = offTargetReadCount;
	}

	@Override
	public String toString() {
		return String.format("%s", getSamRecord().getReadName());
	}

	@Override
	public void analyzeCoverage() {
		// in this analysis, skip and count the reads which are flagged as
		// unmapped and return for getting next read
		if (this.getSamRecord().getReadUnmappedFlag()) {
			readUnmappedCount++;
			return;
		} else {
			// get the offset for the absolute localization
			long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
					(short) (this.getSamRecord().getReferenceIndex() + 1))
					.getOffset();
			long absAlStartOfRead = 0;
			long absAlEndOfRead = 0;
			// Read is a forward Read ...
			if (!this.getSamRecord().getReadNegativeStrandFlag()) {
				absAlStartOfRead = this.getSamRecord().getAlignmentStart()
						+ offset;
				absAlEndOfRead = this.getSamRecord().getAlignmentEnd() + offset;
				// else reverse Read
			} else { // auch bei reversen Reads ist absAlStart < absAlEnd
				absAlStartOfRead = this.getSamRecord().getAlignmentStart()
						+ offset;
				absAlEndOfRead = this.getSamRecord().getAlignmentEnd() + offset;
			}

			/*
			 * create a Read as an IntervalAbs as an inverse Read and use
			 * floor() on the TreeSet to identify the affected IntervalAbs very
			 * fast. We need to use the inverse Read because we want to identify
			 * an IntervalAbs for a Read starting before the IntervalAbs, but
			 * overlaps at their end with the interval.
			 * 
			 * The IntervalAbs is stored in the Reference variable floorInterval
			 * and used in the IntervalAbs-method
			 * readHitsAbsInt(tempRead.getEndAbs(),
			 * tempRead.getStartAbs):boolean with the initial start- and end-
			 * position of the read
			 */
			IntervalAbs tempRead = new IntervalAbs((short) (getSamRecord()
					.getReferenceIndex() + 1), absAlEndOfRead,
					absAlStartOfRead,
					(int) (absAlEndOfRead - absAlStartOfRead), this
							.getSamRecord().getReadName());

			IntervalAbs floorInterval = CheKov.getIntervalTreeSet().floor(
					tempRead);
			if (floorInterval != null) {
				if (floorInterval.readHitsAbsInterval(tempRead.getEndAbs(),
						tempRead.getStartAbs())) {
					// count the reads at least hit one position of the
					// IntervalAbs
					FragmentReadEntry.setOnTargetReadCount(FragmentReadEntry
							.getOnTargetReadCount() + 1);
					// for the moment skip the reverse orientated reads for
					// trouble shooting
					// if (!this.getSamRecord().getReadNegativeStrandFlag())
					this.analyseQuality();
				} else { // count if read is NOT on target, do not analyze
							// quality
					FragmentReadEntry.setOffTargetReadCount(FragmentReadEntry
							.getOffTargetReadCount() + 1);
					// this.analyseQuality();
				}
			}
		}
	}

	@Override
	public void analyseQuality() {
		// CigarString is splitted using StringTokenizer
		String cigarString = getSamRecord().getCigarString();
		StringTokenizer st = new StringTokenizer(cigarString, "MIDNSHPX", true);

		while (st.hasMoreTokens()) {
			cigarTokens.add(st.nextToken());
		}
		int posInRead = 0;
		int posInContig = 0;
		long posInAbsGenom = 0;
		int numberOfBasesAffected = 0;
		int index = 0;
		ReferenceReadPosition rrp = null;
		for (String s : cigarTokens) {
			switch (s) {
			case "H":
				this.setHardClippedBases(this.getHardClippedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				// the corresponding length of the tag is always one field
				// before the tag
				// posInRead = posInRead
				// + Integer.parseInt(cigarTokens.get(index - 1));
				break;
			case "S":
				this.setSoftClippedBases(this.getSoftClippedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));
				// the corresponding length of the tag is always one field
				// before the tag
				// posInRead = posInRead
				// + Integer.parseInt(cigarTokens.get(index - 1));
				break;

			case "D":
				// calculate the total number of by deletion affected bases
				// the corresponding length of the tag is always one field
				// before the tag
				this.setDeletedTaggedBases(this.getDeletedTaggedBases()
						+ Integer.parseInt(cigarTokens.get(index - 1)));

				numberOfBasesAffected = Integer.parseInt(cigarTokens
						.get(index - 1));
				// cumulative position in the read
				posInRead = posInRead + numberOfBasesAffected;
				/*
				 * we want to check, if and what length a putative Homopolymer
				 * is 1. first identify the position in the read we have to the
				 * subtract the length of the affected bases (=
				 * numberOfBasesAffected)
				 */
				posInContig = this.getSamRecord().getAlignmentStart()
						+ posInRead - numberOfBasesAffected;
				// this variable is for identifying the position in the genome
				// using one parameter only, and to create a
				// TargetNucleotidePositionEntry that holds all alterations
				posInAbsGenom = posInContig
						+ ChromosomeOffset.getChromosomeOffsetbyNumber(
								(short) (this.getSamRecord()
										.getReferenceIndex() + 1)).getOffset();
				/*
				 * 2. identify the position in reference genome, Interval AND/OR
				 * IntervalAbs 3. identify and check the missing nucleotide 4.
				 */
				rrp = ReferenceReadPosition.getReferenceReadPositionInstance(
						this.getSamRecord().getReferenceName(), posInContig,
						10, 10);

				/*
				 * create a TargetNucleotidePositionEntry and if new add it to
				 * the ArrayList otherwise update the entry
				 */
				TargetNucleotidePositionEntry tnpe = new TargetNucleotidePositionEntry(
						(short) (this.getSamRecord().getReferenceIndex() + 1),
						posInAbsGenom, AlterationType.Del, rrp.getRefAllel(),
						'-', rrp);
				/*
				 * we need to implement override equals() in
				 * TargetNucleotidePositionEntry
				 */
				// if the altered position does not exist add it to Set
				if (!CheKov.getAlteredNucleotidePositionsEntries().contains(
						tnpe)) {
					TreeSet<TargetNucleotidePositionEntry> tempSet = CheKov
							.getAlteredNucleotidePositionsEntries();
					tempSet.add(tnpe);
					CheKov.setAlteredNucleotidePositionsEntries(tempSet);

					/*
					 * calculate the homoPolymerLength once when
					 * TargetNucleotidePositionEntry is firstly initialized.
					 * this is a stable value for a
					 * TargetNucleotidePositionEntry
					 */
					byte[] tempArray = rrp.getNucleotideInReference();
					int count = 1;
					// count from right
					for (int i = rrp.getIndexOfRefAllele() + 1; i < tempArray.length; i++) {
						if (tempArray[i] == rrp.getRefAllel())
							count++;
						else
							break;
					}
					// count from left
					for (int k = rrp.getIndexOfRefAllele() - 1; k >= 0; k--) {
						if (tempArray[k] == rrp.getRefAllel())
							count++;
						else
							break;
					}
					tnpe.setHomoPolymerLength(count);

				} else {
					/*
					 * update TargetNucleotidePositionEntry (1) find the entry
					 * in TreeSet (2) remove the entry (3) update the entry (4)
					 * add to TreeSet again
					 */
					TargetNucleotidePositionEntry tempEntry = CheKov
							.getAlteredNucleotidePositionsEntries().floor(tnpe);
					// (2) remove the entry
					CheKov.getAlteredNucleotidePositionsEntries().remove(
							tempEntry);
					// (3) Update the entry
					tempEntry.setCoverage(tempEntry.getCoverage() + 1);
					tempEntry.setAltAllelReadCount(tempEntry
							.getAltAllelReadCount() + 1);
					TreeSet<TargetNucleotidePositionEntry> tempSet = CheKov
							.getAlteredNucleotidePositionsEntries();
					// (4) add to TreeSet again
					tempSet.add(tempEntry);
					CheKov.setAlteredNucleotidePositionsEntries(tempSet);
				} // end case "D"
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
		}
		// calculate the rawReadLenghth by adding softClipped and hardClippes
		// bases to effektive Readlength
		this.setRawReadLength(this.getSamRecord().getReadLength()
				+ this.getHardClippedBases() + this.getSoftClippedBases());
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

	public static long getEndTime() {
		return endTime;
	}

	public static void setEndTime(long endTime) {
		FragmentReadEntry.endTime = endTime;
	}

	public static long getStartTime() {
		return startTime;
	}

	public static void setStartTime(long startTime) {
		FragmentReadEntry.startTime = startTime;
	}

}
