import GeneImagePage from "../GeneImagePage";
import defaultImg from "../../assets/imgs/mini_heart_default_plots.png";

export default function ScrnaMiniHeart() {
  return (
    <GeneImagePage
      title="Mini-heart"
      baseUrl="http://128.84.41.80:9029"
      pathPrefix="/genes"
      defaultImageSrc={defaultImg}
    />
  );
}
