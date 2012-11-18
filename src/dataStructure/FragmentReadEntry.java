package dataStructure;

import net.sf.samtools.SAMRecord;

public class FragmentReadEntry extends ReadEntry {

	

	public FragmentReadEntry(SAMRecord samRecord) {
		this.setSamRecord(samRecord);
	}

	public SAMRecord getSamRecord() {
		return super.getSamRecord();
	}

	public void setSamRecord(SAMRecord samRecord) {
		super.setSamRecord(samRecord);
	}

	@Override
	public String toString() {
		return String.format("%s",getSamRecord().getReadName());
	}

	@Override
	public void analyzeCoverage() {
		// TODO Auto-generated method stub
		
	}
	
	


}
