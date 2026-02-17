import { create } from 'zustand';

type RightTab = 'variables' | 'arrays' | 'config' | 'debug' | 'diagrams';
type BottomTab = 'console' | 'sensitivity' | 'diagrams';

interface UIState {
  theme: 'light' | 'dark';
  rightTab: RightTab;
  bottomTab: BottomTab;
  bottomPanelOpen: boolean;
  splitPosition: number; // percentage for left pane

  setTheme: (theme: 'light' | 'dark') => void;
  toggleTheme: () => void;
  setRightTab: (tab: RightTab) => void;
  setBottomTab: (tab: BottomTab) => void;
  setBottomPanelOpen: (open: boolean) => void;
  toggleBottomPanel: () => void;
  setSplitPosition: (pos: number) => void;
}

export const useUIStore = create<UIState>((set) => ({
  theme: (localStorage.getItem('coolsolve-theme') as 'light' | 'dark') || 'light',
  rightTab: 'variables',
  bottomTab: 'console',
  bottomPanelOpen: true,
  splitPosition: 50,

  setTheme: (theme) => {
    localStorage.setItem('coolsolve-theme', theme);
    set({ theme });
  },
  toggleTheme: () =>
    set((s) => {
      const next = s.theme === 'light' ? 'dark' : 'light';
      localStorage.setItem('coolsolve-theme', next);
      return { theme: next };
    }),
  setRightTab: (tab) => set({ rightTab: tab }),
  setBottomTab: (tab) => set({ bottomTab: tab, bottomPanelOpen: true }),
  setBottomPanelOpen: (open) => set({ bottomPanelOpen: open }),
  toggleBottomPanel: () => set((s) => ({ bottomPanelOpen: !s.bottomPanelOpen })),
  setSplitPosition: (pos) => set({ splitPosition: pos }),
}));
