package algorithm;

public enum AppParameterCLI {
	BED("b", "bed", true, "read the bed-File"),
	BAM("a", "bam", true, "read the bam-File"),
	MISSED("m", "missed", true, "output-File to write the missed intervals"),
	OUT("o","out", true, "output-File to write the coverages of intervals"),
	REF("r", "ref", true, "read the reference-File"),
	HOM("h", "hom", true, "output-File to write the homopolymers"),
	THRESHOLD("t", "threshold", true, "the threshold Integer plus/minus region of the interval"),
	UNDERREPRESENTED("u","under",true, "the cutoff coverage for identifying an underrepresented"),
	UNDERFILE("uf", "underFile",true, "the filename for the underrepresented")
	;
	final String shortName;
	final String name;
	final boolean hasArgs;
	final String helpText;

	private AppParameterCLI(final String shortName, final String name,
			final boolean hasArgs, final String helpText) {
		this.shortName = shortName;
		this.name = name;
		this.hasArgs = hasArgs;
		this.helpText = helpText;
		
	}
}
