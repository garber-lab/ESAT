package broad.core.hmm;

import java.util.ArrayList;
import java.util.List;

public class BadModelException extends Exception {

	private static final long serialVersionUID = 712946940195037285L;
	List<String> problems;
	
	public BadModelException(MarkovModel m) {
		if(m == null) {
			problems = new ArrayList<String>(1);
			problems.add("MarkovModel is null!");
		} else {
			problems = m.getConsistencyProblems();
		}
	}

	public String getMessage() {
		StringBuilder sb = new StringBuilder();
		if(problems.size() == 0) {
			sb.append("Unknown problem");
		} else {
			for(int i = 0; i < problems.size(); i++) {
				sb.append("\n").append(i).append(". ").append(problems.get(i));
			}
		}
		
		return sb.toString();
	}
}
