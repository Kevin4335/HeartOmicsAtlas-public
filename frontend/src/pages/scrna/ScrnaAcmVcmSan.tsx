import GeneImagePage from "../GeneImagePage";
import defaultImg from "../../assets/imgs/ACM_VCM_SAN_default_plots.png";

export default function ScrnaAcmVcmSan() {
  return (
    <GeneImagePage
      title="ACM_VCM_SAN"
      baseUrl="http://128.84.41.80:9027"
      pathPrefix="/genes"
      defaultImageSrc={defaultImg}
    />
  );
}
