import GeneImagePage from "./GeneImagePage";
import stDefault from "../assets/imgs/Multiomics_default_plots.png";

export default function SpatialTranscriptomics() {
  return (
    <GeneImagePage
      title="scMultiomics"
      baseUrl="http://128.84.41.80:9026"
      pathPrefix="/genes"
      defaultImageSrc={stDefault}
    />
  );
}
