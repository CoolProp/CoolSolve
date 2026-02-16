import { useRef, useCallback, useEffect, useState } from 'react';

interface SplitPaneProps {
  direction: 'horizontal' | 'vertical';
  defaultSize: number;       // initial size in % for the first pane
  minSize?: number;           // minimum size in px
  maxSize?: number;           // maximum size in px
  children: [React.ReactNode, React.ReactNode];
  className?: string;
  onResize?: (size: number) => void;
}

export default function SplitPane({
  direction,
  defaultSize,
  minSize = 100,
  maxSize,
  children,
  className = '',
  onResize,
}: SplitPaneProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState(defaultSize);
  const dragging = useRef(false);

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    dragging.current = true;
    document.body.style.cursor = direction === 'horizontal' ? 'col-resize' : 'row-resize';
    document.body.style.userSelect = 'none';
  }, [direction]);

  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      if (!dragging.current || !containerRef.current) return;

      const rect = containerRef.current.getBoundingClientRect();
      let newSize: number;

      if (direction === 'horizontal') {
        const offset = e.clientX - rect.left;
        newSize = (offset / rect.width) * 100;
      } else {
        const offset = e.clientY - rect.top;
        newSize = (offset / rect.height) * 100;
      }

      // Apply min/max constraints
      const containerPx = direction === 'horizontal' ? rect.width : rect.height;
      const minPct = (minSize / containerPx) * 100;
      const maxPct = maxSize ? (maxSize / containerPx) * 100 : 90;

      newSize = Math.max(minPct, Math.min(maxPct, newSize));
      setSize(newSize);
      onResize?.(newSize);
    };

    const handleMouseUp = () => {
      if (dragging.current) {
        dragging.current = false;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
      }
    };

    window.addEventListener('mousemove', handleMouseMove);
    window.addEventListener('mouseup', handleMouseUp);
    return () => {
      window.removeEventListener('mousemove', handleMouseMove);
      window.removeEventListener('mouseup', handleMouseUp);
    };
  }, [direction, minSize, maxSize, onResize]);

  const isHorizontal = direction === 'horizontal';

  return (
    <div
      ref={containerRef}
      className={`split-pane ${isHorizontal ? 'split-horizontal' : 'split-vertical'} ${className}`}
      style={{
        display: 'flex',
        flexDirection: isHorizontal ? 'row' : 'column',
        flex: 1,
        minHeight: 0,
        minWidth: 0,
        overflow: 'hidden',
      }}
    >
      {/* First pane */}
      <div
        className="split-pane-first"
        style={{
          [isHorizontal ? 'width' : 'height']: `${size}%`,
          minWidth: isHorizontal ? minSize : undefined,
          minHeight: !isHorizontal ? minSize : undefined,
          overflow: 'hidden',
          display: 'flex',
          flexDirection: 'column',
        }}
      >
        {children[0]}
      </div>

      {/* Divider */}
      <div
        className={`split-divider ${isHorizontal ? 'split-divider-h' : 'split-divider-v'}`}
        onMouseDown={handleMouseDown}
      />

      {/* Second pane */}
      <div
        className="split-pane-second"
        style={{
          flex: 1,
          minWidth: isHorizontal ? minSize : undefined,
          minHeight: !isHorizontal ? minSize : undefined,
          overflow: 'hidden',
          display: 'flex',
          flexDirection: 'column',
        }}
      >
        {children[1]}
      </div>
    </div>
  );
}
