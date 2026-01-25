import GeneImagePage from "./GeneImagePage";
import stDefault from "../assets/imgs/Spatial_default_plots.png";

export default function SpatialTranscriptomics() {
  return (
    <GeneImagePage
      title="Spatial Transcriptomics"
      baseUrl="http://128.84.41.80:9025"
      pathPrefix="/genes"
      defaultImageSrc={stDefault}
    />
  );
}
