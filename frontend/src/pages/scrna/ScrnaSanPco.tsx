import GeneImagePage from "../GeneImagePage";
import defaultImg from "../../assets/imgs/SAN_PCO_default_plots.png";

export default function ScrnaSanPco() {
  return (
    <GeneImagePage
      title="SAN-PCO"
      baseUrl="http://128.84.41.80:9028"
      pathPrefix="/genes"
      defaultImageSrc={defaultImg}
    />
  );
}
