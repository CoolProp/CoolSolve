import { useState, useRef, useCallback, type ReactNode } from 'react';

interface TooltipProps {
  /** Content rendered inside the tooltip bubble */
  content: ReactNode;
  /** The trigger element (e.g. an info icon) */
  children: ReactNode;
}

/**
 * A viewport-aware tooltip that appears above its trigger element.
 * Automatically clamps horizontal position so it never overflows
 * the left or right edge of the viewport.
 */
export default function Tooltip({ content, children }: TooltipProps) {
  const [visible, setVisible] = useState(false);
  const [pos, setPos] = useState({ left: 0, bottom: 0 });
  const wrapRef = useRef<HTMLSpanElement>(null);

  const show = useCallback(() => {
    const el = wrapRef.current;
    if (!el) return;
    const rect = el.getBoundingClientRect();
    const tooltipW = 280;
    const gap = 6;
    const margin = 4;                       // min distance from viewport edge

    // Ideal: centred on the icon
    let left = rect.left + rect.width / 2 - tooltipW / 2;
    // Clamp so the tooltip stays within the viewport
    left = Math.max(margin, Math.min(left, window.innerWidth - tooltipW - margin));

    setPos({ left, bottom: window.innerHeight - rect.top + gap });
    setVisible(true);
  }, []);

  const hide = useCallback(() => setVisible(false), []);

  return (
    <span
      ref={wrapRef}
      className="config-tooltip-wrap"
      onMouseEnter={show}
      onMouseLeave={hide}
    >
      {children}
      {visible && (
        <span
          className="config-tooltip"
          style={{ left: pos.left, bottom: pos.bottom }}
        >
          {content}
        </span>
      )}
    </span>
  );
}
