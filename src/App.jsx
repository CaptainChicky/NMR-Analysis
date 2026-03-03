import { useState } from "react";
import CarbonNMR from "./carbonNMR";
import ProtonNMR from "./protonNMR";

export default function App() {
	const [nucleus, setNucleus] = useState("13C");

	return nucleus === "1H" ? (
		<div>
			<div style={{ padding: "12px 16px 0", background: "#0f172a", display: "flex", gap: 8, alignItems: "center" }}>
				<button onClick={() => setNucleus("13C")} style={{ padding: "6px 16px", borderRadius: 8, border: "none", cursor: "pointer", background: "#1e293b", color: "#94a3b8", fontWeight: 600, fontSize: 13 }}>¹³C</button>
				<button style={{ padding: "6px 16px", borderRadius: 8, border: "none", background: "#3b82f6", color: "#fff", fontWeight: 600, fontSize: 13 }}>¹H</button>
			</div>
			<ProtonNMR />
		</div>
	) : (
		<div>
			<div style={{ padding: "12px 16px 0", background: "#0f172a", display: "flex", gap: 8, alignItems: "center" }}>
				<button style={{ padding: "6px 16px", borderRadius: 8, border: "none", background: "#3b82f6", color: "#fff", fontWeight: 600, fontSize: 13 }}>¹³C</button>
				<button onClick={() => setNucleus("1H")} style={{ padding: "6px 16px", borderRadius: 8, border: "none", cursor: "pointer", background: "#1e293b", color: "#94a3b8", fontWeight: 600, fontSize: 13 }}>¹H</button>
			</div>
			<CarbonNMR />
		</div>
	);
}