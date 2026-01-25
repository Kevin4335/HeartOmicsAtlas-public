declare module "react-zoom-pan-pinch" {
  import type { CSSProperties, ReactNode } from "react";

  export interface ReactZoomPanPinchContent {
    resetTransform: () => void;
    zoomIn: () => void;
    zoomOut: () => void;
  }

  export interface TransformWrapperProps {
    children: (utils: ReactZoomPanPinchContent) => ReactNode;
    initialScale?: number;
    minScale?: number;
    maxScale?: number;
    centerOnInit?: boolean;
    limitToBounds?: boolean;
    wheel?: { step?: number };
    doubleClick?: { disabled?: boolean };
    panning?: { velocityDisabled?: boolean };
    [key: string]: unknown;
  }

  export interface TransformComponentProps {
    wrapperStyle?: CSSProperties;
    contentStyle?: CSSProperties;
    children?: ReactNode;
  }

  export function TransformWrapper(props: TransformWrapperProps): JSX.Element;
  export function TransformComponent(props: TransformComponentProps): JSX.Element;
}
