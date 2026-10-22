const Highlight = ({children, color}) => (
  <span
    style={{
      backgroundColor: color || 'var(--highlight-background-color)',
      borderRadius: '2px',
      color: 'var(--highlight-text-color)',
      padding: '0.2rem',
    }}>
    {children}
  </span>
);

const HighlightOptionalArg = ({children}) => (
  <span
    style={{
      backgroundColor: 'var(--highlight-optional-arg-background-color)',
      borderRadius: '2px',
      color: 'var(--highlight-optional-arg-text-color)',
      padding: '0.2rem',
    }}>
    {children}
  </span>
);

export { Highlight, HighlightOptionalArg };
